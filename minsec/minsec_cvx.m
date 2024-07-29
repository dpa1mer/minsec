function [S, C, f, field] = minsec_cvx(meshData, N, lambda, degree, fiberRadius)

nv = meshData.nv;
ne = meshData.ne;
nf = meshData.nf;
nb = meshData.nb;
bdryIdx = meshData.bdryIdx;

bdryEdgeIdx = meshData.bdryEdgeIdx;

%% Curvature
scaledCurvature = degree * meshData.angleDefect / (2 * pi);
if nb > 0
    scaledCurvature(meshData.bdryIdx) = 0;
end
vtxPointwiseScaledCurvature = meshData.star0lump \ scaledCurvature;

% Redistribute curvature to edge centers as we're using the non-conforming
% Crouzeix-Raviart finite elements for our singularities
pointwiseScaledCurvature = (vtxPointwiseScaledCurvature(meshData.edges(:, 1)) + ...
                            vtxPointwiseScaledCurvature(meshData.edges(:, 2))) / 2;

%% Operators
[k, ik, ik_pinv, d_vert] = fourier_ops(N, false, true);
K = length(k);

FT = fft(eye(N), [], 1);
FT = FT(1:K, :);
FT_blocks = kron(FT, speye(9 * nf));

% Note that Fourier components get parallel transported with an *inverse*
% phase compared to vectors
[G_intr, M_intr, G_nc, M_nc] = intrinsic_ops(meshData);
[V2F, ~, ~, F2C] = covariant_ops(meshData, -degree * k);

AreaFx2 = kron(spdiags(meshData.areas, 0, nf, nf), speye(2));

Horiz = kron(speye(3 * nf), sparse([1 0 0; 0 1 0]));
Vert = kron(speye(3 * nf), sparse([0 0 1]));
assert(norm(Horiz' * Horiz + Vert' * Vert - speye(3 * 3 * nf), 'fro') == 0);

G_total = (kron(speye(K), Horiz.' * kron(F2C, speye(2)) * G_intr) + ...
           kron(spdiags(d_vert, 0, K, K), Vert.')) * V2F;

J = kron(speye(nf), [0 -1; 1 0]); % rotation by pi/2 - acts like Hodge star on vectors
JG_nc = J * G_nc;
F2Cx2 = kron(F2C, [1 0; 0 1; 0 0]);
L_nc = JG_nc' * AreaFx2 * JG_nc;

%% Boundary
% Break symmetry by fixing the values of f on the boundary
if nb > 0
    [bdryAngles, ~, BdryInt, bdryIntVal] = field_bdry_conds(meshData, degree);

    f_hifreq_fixed = ((N-1) / (2 * pi)) * (ik_pinv .* (1 - abs(k) ./ (N / 2))).';
    f_hifreq_fixed = f_hifreq_fixed(2:end);
    f_hifreq_fixed = f_hifreq_fixed .* exp(-ik(2:end).' .* bdryAngles);
end

%% Define optimization problem
tau_bar = repmat([0; 0; 1 / (2*pi)], 1, 3, nf);

g = spdiags([fiberRadius; fiberRadius; 1], 0, 3, 3);

cvx_begin
    variable S(3, 3, nf, N)
    variable C(ne, 1)
    variable f_ft(nv, K) complex
    variable phi(ne, 1)

    df_ft = G_total * f_ft(:);
    S_ft = reshape(df_ft, 3, 3, nf, K);
    S_ft(:, :, :, 1) = S_ft(:, :, :, 1) + reshape(N * (F2Cx2 * JG_nc * phi + tau_bar(:)), 3, 3, nf);

    minimize ((2 * pi / N) * sum(sum(M_intr * reshape(norms(g * reshape(S, 3, []), 2, 1), 3 * nf, N), 1), 2) + ...
             lambda * sum(M_nc * abs(C)))
    subject to
        FT_blocks * S(:) == S_ft(:);
        L_nc * phi == M_nc * (C - pointwiseScaledCurvature);
        S(3, :, :, :) >= 0;

        phi(bdryEdgeIdx) == 0;
        BdryInt * reshape(sum(S_ft(1:2, :, :, 1), 2) / 3, [], 1) == N * bdryIntVal;
        f_ft(bdryIdx, 2:end) == f_hifreq_fixed;
cvx_end

f = ifft(f_ft, [], 2, 'symmetric');
field = -1i * f_ft(:, end);

end