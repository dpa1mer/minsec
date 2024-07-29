function [S, C, f, field, admmData] = minsec(meshData, N, lambda, degree, fiberRadius, options, vizOpts)

arguments
    meshData struct
    N           (1, 1) double {mustBeInteger, mustBePositive} = 64
    lambda      (:, 1) double {mustBePositive}                = 1
    degree      (1, 1) double {mustBeInteger, mustBePositive} = 1
    fiberRadius (1, 1) double {mustBePositive}                = 1

    options.MaxIters                     (:, 1) double {mustBeInteger,mustBePositive} = 10000
    options.UseGPU                       (1, 1) logical = (gpuDeviceCount > 0)
    options.UseDiscreteVerticalOperators (1, 1) logical = false
    options.FixedSingOnly                (1, 1) logical = false
    options.FixedSingularities           (:, 1) double {mustBeInteger} = []
    options.ExcludedSingularities        (:, 1) double {mustBeInteger} = []

    vizOpts.Visualize (1, 1) logical                               = false
    vizOpts.Scale     (1, 1) double {mustBePositive}               = 1
    vizOpts.Repeat    (1, 1) double {mustBePositive,mustBeInteger} = 1
    vizOpts.Frequency (1, 1) double {mustBePositive,mustBeInteger} = 100
    vizOpts.Shift     (1, 1) double                                = 1
    vizOpts.SubsampleRate (1, 1) double { ...
        mustBePositive, ...
        mustBeLessThanOrEqual(vizOpts.SubsampleRate, 1)} = 1
    vizOpts.Directory (1, 1) string
end

if options.UseGPU
    data_type = "gpuArray";
else
    data_type = "double";
end


%% Convenience definitions

nv = meshData.nv;
nf = meshData.nf;
nb = meshData.nb;
ni = nv - nb;
ne = meshData.ne;
intEdgeIdx = meshData.intEdgeIdx;
nie = length(intEdgeIdx);

fiberHeight = 2 * pi * fiberRadius;

% Handle spatially varying regularization parameter
assert(size(lambda, 1) == ne || size(lambda, 1) == 1);
if size(lambda, 1) == ne
    lambda = lambda(intEdgeIdx);
end


%% Global operators

[k, ik, ik_pinv, d_vert] = fourier_ops(N, options.UseDiscreteVerticalOperators, true);
% K: Number of frequencies. Since the current is real, we only need half as
% many frequencies as discrete points due to conjugate symmetry.
K = length(k);

% Note that Fourier components get parallel transported with an *inverse*
% phase compared to vectors
[G_intr, M_intr, G_nc, M_nc] = intrinsic_ops(meshData);
[V2F, ~, ~, F2C] = covariant_ops(meshData, -degree * k);
AreaFx2 = kron(spdiags(meshData.areas, 0, nf, nf), speye(2));

% Pointwise inner product for 1-forms on the bundle
gij_bundle = spdiags([1; 1; fiberRadius^(-2)], 0, 3, 3);
gij_sqrt_bundle = [1; 1; fiberRadius^(-1)];

% L^2 for 1-forms on the bundle
M_bundle = (fiberHeight / N) * kron(M_intr, gij_bundle);

Horiz = kron(speye(3 * nf), sparse([1 0 0; 0 1 0]));
Vert = kron(speye(3 * nf), sparse([0 0 1]));
assert(norm(Horiz' * Horiz + Vert' * Vert - speye(3 * 3 * nf), 'fro') == 0);

G_total = (kron(speye(K), Horiz.' * kron(F2C, speye(2)) * G_intr) + ...
           kron(spdiags(d_vert, 0, K, K), Vert.')) * V2F;
G_horiz = Horiz.' * kron(F2C, speye(2)) * G_intr;
L_total = G_total' * (kron(speye(K), M_bundle) * G_total);


%% Boundary conditions
if nb > 0
    [bdryAngles, ~, BdryInt, bdryIntVal] = field_bdry_conds(meshData, degree);
end


%% High-frequency part
f_hifreq_fixed = ((N - 1) / (2 * pi)) * (ik_pinv .* (1 - abs(k) ./ (N / 2))).';
f_hifreq_fixed = f_hifreq_fixed(2:end);
if nb > 0
    FixedPart = sparse(1:nb, meshData.bdryIdx, 1, nb, nv);
    VarPart = sparse(1:ni, meshData.intIdx, 1, ni, nv);

    f_hifreq_fixed = f_hifreq_fixed .* exp(-ik(2:end).' .* bdryAngles);
else
    FixedPart = sparse(1, 1, 1, 1, nv);
    VarPart = sparse(1:(nv - 1), 2:nv, 1, nv - 1, nv);
end
assert(norm(FixedPart' * FixedPart + VarPart' * VarPart - speye(nv), 'fro') == 0);

L_hifreq = L_total((nv+1):end, (nv+1):end);
f_hifreq_fixed = FixedPart' * f_hifreq_fixed;
f_hifreq_fixed_rhs = VarPart * reshape(L_hifreq * f_hifreq_fixed(:), nv, K - 1);

f_hifreq_var = [];

% The main system matrix for the high-frequency component
HiFreqVarPart = kron(speye(K - 1), VarPart);
L_hifreq = HiFreqVarPart * L_hifreq * HiFreqVarPart';
L_hifreq = (L_hifreq + L_hifreq') / 2;

hifreq_solver = schur_solver(L_hifreq, [], options.UseGPU);


%% Low-frequency part
% Non-conforming finite elements for Biot-Savart field computation

% Fix f0 value at first vertex to zero to kill kernel
G0 = G_intr * V2F(1:(3*nf), 2:nv);
L0 = G0' * AreaFx2 * G0;

J = kron(speye(nf), [0 -1; 1 0]); % rotation by pi/2 - acts like Hodge star on vectors
JG_tngt = J * G_nc(:, intEdgeIdx);
L_tngt = JG_tngt' * AreaFx2 * JG_tngt;
M_tngt = M_nc(intEdgeIdx, intEdgeIdx);
BiL_tngt = L_tngt' * (M_tngt \ L_tngt);

if nb > 0
    Constraint_lofreq = BdryInt * [G0, JG_tngt];
    constraintVal = bdryIntVal;
else
    % Set last edge value of g0 to zero [0, ..., 0, 1]
    Constraint_lofreq = sparse(1, nv - 1 + ne, 1, 1, nv - 1 + ne);
    constraintVal = 0;
end

mu = 1;
nu = 1;
[lofreq_solver, mu_scaled] = recalculateFreq0Operators(mu, nu);


%% Gauss curvature, expressed in Crouzeix-Raviart basis
scaledCurvature = degree * meshData.angleDefect / (2 * pi);
if nb > 0
    scaledCurvature(meshData.bdryIdx) = 0;
end
vtxPointwiseScaledCurvature = meshData.star0lump \ scaledCurvature;

% Redistribute curvature to edge centers as we're using the non-conforming
% Crouzeix-Raviart finite elements for our singularities
pointwiseScaledCurvature = (vtxPointwiseScaledCurvature(meshData.edges(:, 1)) + ...
                            vtxPointwiseScaledCurvature(meshData.edges(:, 2))) / 2;
pointwiseScaledCurvature = pointwiseScaledCurvature(intEdgeIdx);


%% Singularity constraints
CFixedIdx = abs(options.FixedSingularities);
[~, CFixedIdx] = ismember(CFixedIdx, intEdgeIdx);
assert(all(CFixedIdx > 0));
CFixedVal = M_tngt(CFixedIdx, CFixedIdx) \ sign(options.FixedSingularities);

% Handle singularity exclusion zones
[~, CZeroIdx] = ismember(options.ExcludedSingularities, intEdgeIdx);
CZeroIdx = nonzeros(CZeroIdx);

if options.FixedSingOnly
    CZeroIdx = setdiff(1:nie, CFixedIdx);
end


%% Initialize ADMM
tau_bar = repmat([0; 0; 1 / (2*pi)], 1, 3, nf, N);

S = zeros(3, 3, nf, N, data_type);
Sbar = S;

Cbar = zeros(nie, 1, data_type);
C = Cbar;
C(CFixedIdx) = CFixedVal;
Cbar(CFixedIdx) = CFixedVal;

w = zeros(size(tau_bar), data_type);
z = zeros(size(C), data_type);


%% Move stuff to GPU if necessary
if options.UseGPU
    d_vert = gpuArray(d_vert);
    tau_bar = gpuArray(tau_bar);
    pointwiseScaledCurvature = gpuArray(pointwiseScaledCurvature);
    constraintVal = gpuArray(constraintVal);
    V2F = gpuArray(V2F);
    Vert = gpuArray(Vert);
    G_horiz = gpuArray(G_horiz);
    VarPart = gpuArray(VarPart);
    f_hifreq_fixed = gpuArray(f_hifreq_fixed);
    f_hifreq_fixed_rhs = gpuArray(f_hifreq_fixed_rhs);
    M_intr = gpuArray(M_intr);
    M_bundle = gpuArray(M_bundle);
    gij_sqrt_bundle = gpuArray(gij_sqrt_bundle);
    G0 = gpuArray(G0);
    AreaFx2 = gpuArray(AreaFx2);
    L0 = gpuArray(L0);
    JG_tngt = gpuArray(JG_tngt);
    L_tngt = gpuArray(L_tngt);
    M_tngt = gpuArray(M_tngt);
    BiL_tngt = gpuArray(BiL_tngt);
    lambda = gpuArray(lambda);
    CFixedIdx = gpuArray(CFixedIdx);
    CFixedVal = gpuArray(CFixedVal);
    CZeroIdx = gpuArray(CZeroIdx);
    shrink_i = gpuArray.colon(1, 3).';
    shrink_j = gpuArray.colon(1, 3 * nf * N);
end


%% Optional visualization setup

if vizOpts.Visualize
    fig = figure(Color='white'); hold on;
    fig.Position(3:4) = [1024 1024];
    ptc1 = plot_section_volume(meshData, gij_sqrt_bundle .* S, Scale=vizOpts.Scale, ...
                               Repeat=vizOpts.Repeat, Shift=vizOpts.Shift);
    ptc2 = plot_section(meshData, gij_sqrt_bundle .* S, Degree=degree, SubsampleRate=vizOpts.SubsampleRate);
    view(3);
    axis image vis3d off;
    cameratoolbar;

    saveViz = false;
    if isfield(vizOpts, 'Directory')
        saveViz = true;
        mkdir(vizOpts.Directory);
    end
end


%% ADMM

admmData = struct;

SOld = 0; COld = 0;
tic;
for iter = 1:options.MaxIters
    [Sbar, Cbar, f, f_hifreq_var] = global_step(S, C, w, z, f_hifreq_var);
    [S, C] = local_step(Sbar, Cbar, w, z);
    [w, z, pRes1, pRes2] = dual_update(w, z, S, Sbar, C, Cbar);

    pRes1 = sqrt(pRes1(:)' * reshape(M_bundle * reshape(pRes1, 9 * nf, N), [], 1));
    pRes2 = sqrt(pRes2(:)' * M_tngt * pRes2(:));
    fprintf('Primal Residual:\t||S - Sbar - tau_bar|| = %d\t||C - Cbar|| = %d\n', ...
            pRes1, pRes2);

    dRes1 = SOld - S;
    dRes1 = mu * sqrt(dRes1(:)' * reshape(M_bundle * reshape(dRes1, 9 * nf, N), [], 1));

    dRes2 = COld - C;
    dRes2 = nu * sqrt(dRes2(:)' * M_tngt * dRes2(:));
    fprintf('Dual Residual:\t\t||SOld - S|| = %d\t\t||COld - C|| = %d\n', ...
            dRes1, dRes2);

    objVal = ...
        sum(M_tngt * (lambda .* abs(C))) + ...
        sum(M_intr * reshape(vecnorm(gij_sqrt_bundle .* S, 2, 1), 3 * nf, N), 'all') * (fiberHeight / N);
    fprintf('Objective value:\t%d\n', objVal);
    fprintf('Mu:\t%d\tNu:\t%d\n', mu, nu);

    admmData(iter).iter = iter;
    admmData(iter).pRes1 = pRes1;
    admmData(iter).pRes2 = pRes2;
    admmData(iter).dRes1 = dRes1;
    admmData(iter).dRes2 = dRes2;
    admmData(iter).objVal = objVal;
    admmData(iter).mu = mu;
    admmData(iter).nu = nu;
    admmData(iter).wallTime = toc;

    if vizOpts.Visualize && mod(iter, vizOpts.Frequency) == 1
        % Update plot
        plot_section_volume(meshData, gij_sqrt_bundle .* S, ptc1, Scale=vizOpts.Scale, ...
                            Repeat=vizOpts.Repeat, Shift=vizOpts.Shift);
        plot_section(meshData, gij_sqrt_bundle .* S, ptc2, Degree=degree, SubsampleRate=vizOpts.SubsampleRate);
        drawnow;
        if saveViz
            filename = fullfile(vizOpts.Directory, sprintf('volume_%.5d.png', iter));
            exportgraphics(fig, filename, Resolution=150);
        end
    end

    % Convergence criterion
    if pRes1 < 5e-4 && dRes1 < 5e-4 && pRes2 < 5e-4 && dRes2 < 5e-4
        break;
    end

    if mod(iter, 10) == 0
        % Adaptively update the ADMM regularization parameters
        % Also update the scaled dual variables accordingly
        updateMuNu = false;
        if pRes1 > 10 * dRes1
            mu = 2 * mu;
            w = 0.5 * w;
            updateMuNu = true;
        elseif dRes1 > 10 * pRes1
            mu = 0.5 * mu;
            w = 2 * w;
            updateMuNu = true;
        end
    
        if ~options.FixedSingOnly
            if pRes2 > 10 * dRes2 && dRes2 > 0
                nu = 2 * nu;
                z = 0.5 * z;
                updateMuNu = true;
            elseif dRes2 > 10 * pRes2
                nu = 0.5 * nu;
                z = 2 * z;
                updateMuNu = true;
            end
        end
    
        if updateMuNu
            [lofreq_solver, mu_scaled] = recalculateFreq0Operators(mu, nu);
        end
    end

    SOld = S; COld = C;
end


%% Extract field
f_ft = fft(f, [], 2);
field = -1i * f_ft(:, end);

if vizOpts.Visualize && saveViz
    writetable(struct2table(admmData), fullfile(vizOpts.Directory, 'admm_data.csv'));
end


%% Subroutines
function [Sbar, Cbar, f, f_hifreq_var] = global_step(S, C, w, z, f_hifreq_var)
    upstairs = fft(w + S - tau_bar, [], 4);
    % We only need to keep track of one conjugate half. Also set the N/2
    % frequency to zero
    upstairs = upstairs(:, :, :, 1:K);

    downstairs = z + C - pointwiseScaledCurvature;

    upstairs_rescaled = M_bundle * reshape(upstairs, 9 * nf, K);
    div_upstairs = reshape(V2F' * reshape(Vert * (upstairs_rescaled .* d_vert') + ...
                                          G_horiz' * upstairs_rescaled, [], 1), nv, K);

    % Solve for DC component
    [Cbar, f0, g0] = global_lofreq(upstairs, downstairs);

    % Solve for higher-frequency components of f
    [f_hifreq, f_hifreq_var] = global_hifreq(div_upstairs, f_hifreq_var);

    f = [0; f0; f_hifreq(:)];
    f_transp = reshape(V2F * f, 3 * nf, K);
    df = Vert.' * (f_transp .* d_vert.') + G_horiz * f_transp;

    Sbar = reshape(df, 3, 3, nf, K);
    Sbar(1:2, :, :, 1) = Sbar(1:2, :, :, 1) + reshape(JG_tngt * g0, 2, 1, nf, 1);

    % Reconstitute the spatial domain values using conjugate symmetry.
    Sbar = ifft(Sbar, N, 4, 'symmetric');
    f = ifft(reshape(f, nv, K), N, 2, 'symmetric');
end

function [f_hifreq, f_hifreq_var] = global_hifreq(div_upstairs, f_hifreq_var)
    f_rhs = VarPart * div_upstairs(:, 2:end) - f_hifreq_fixed_rhs;
    f_hifreq_var = reshape(hifreq_solver.solve(f_rhs(:), f_hifreq_var), [], K - 1);
    f_hifreq = VarPart' * f_hifreq_var + f_hifreq_fixed;
end

function [Cbar, f0, g0] = global_lofreq(upstairs, downstairs)
    projHoriz = (1 / (3 * N)) * reshape(sum(upstairs(1:2, :, :, 1), 2), 2 * nf, 1); % 0-frequency component
    f0_rhs = (mu_scaled * G0' * AreaFx2 * projHoriz);
    g0_rhs = (mu_scaled * JG_tngt' * AreaFx2 * projHoriz + nu * L_tngt * downstairs);

    % Solve Schur complement system to enforce boundary conditions
    fg = lofreq_solver.solve([f0_rhs; g0_rhs], constraintVal);
    f0 = fg(1:(nv-1));
    g0 = fg(nv:end);

    Cbar = M_tngt \ (L_tngt * g0) + pointwiseScaledCurvature;

    f0 = N * f0;
    g0 = N * g0;
end

function [S, C] = local_step(Sbar, Cbar, w, z)
    S = tau_bar + Sbar - w;

    if options.UseGPU
        S = reshape(S, 3, 3 * nf * N);
        S = arrayfun(@shrink_S_gpu, shrink_i, shrink_j);
        S = reshape(S, 3, 3, nf, N);
    else
        S(3, :, :, :) = max(0, S(3, :, :, :)); % vertical component should be positive
        S = S .* max(0, 1 - (mu * vecnorm(gij_sqrt_bundle .* S, 2, 1)).^(-1));
    end
    
    C = Cbar - z;
    C = C .* max(0, 1 - lambda ./ (nu .* abs(C)));

    % Fixed singularities
    C(CFixedIdx) = CFixedVal;

    % Excluded singularities
    C(CZeroIdx) = 0;

    function s = shrink_S_gpu(i, j)
        if i == 3
            s_val = max(0, S(i, j));
        else
            s_val = S(i, j);
        end
        norm_s = sqrt((gij_sqrt_bundle(1) * S(1, j))^2 + ...
                      (gij_sqrt_bundle(2) * S(2, j))^2 + ...
                      (gij_sqrt_bundle(3) * max(0, S(3, j)))^2);
        s = s_val * max(0, 1 - (mu * norm_s)^(-1));
    end
end

function [w, z, pRes1, pRes2] = dual_update(w, z, S, Sbar, C, Cbar)
    pRes1 = S - tau_bar - Sbar;
    w = w + pRes1;

    pRes2 = C - Cbar;
    z = z + pRes2;
end

function [lofreq_solver, mu_scaled] = recalculateFreq0Operators(mu, nu)
    mu_scaled = fiberHeight * mu;

    Op_f0 = gather(mu_scaled * L0);
    Op_g0 = gather(mu_scaled * L_tngt + nu * BiL_tngt);
    System_lofreq = blkdiag(Op_f0, Op_g0);

    lofreq_solver = schur_solver(System_lofreq, Constraint_lofreq, options.UseGPU);
end

end