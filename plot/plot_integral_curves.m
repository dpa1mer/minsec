% plot_integral_curves(meshData, z, degree, ncurves)
% Plot integral curves of a directional field. Assumes a flat domain.
function hdl = plot_integral_curves(meshData, z, degree, ncurves)

if nargin < 4
    ncurves = 2000;
end

nv = meshData.nv;
ne = meshData.ne;
nf = meshData.nf;
verts = meshData.verts;
faces = meshData.faces;

G_intr = intrinsic_ops(meshData);
V2F = covariant_ops(meshData, degree);

face2face = spones(meshData.d1(:, meshData.face2edge)) - [speye(nf) speye(nf) speye(nf)];
[mask, face2face] = max(face2face, [], 1);
face2face(~mask) = nan;
face2face = reshape(face2face, nf, 3);

f2eTransport = sparse(meshData.face2edge, repmat((1:nf).', 1, 3), meshData.f2eTransport, ne, nf);
f2fTransport = f2eTransport' * f2eTransport;
f2fTransport = f2fTransport - spdiags(diag(f2fTransport), 0, nf, nf) + speye(nf, nf);

if size(z, 1) == nv
    % parallel-transport to face corners
    z = reshape(V2F * z, 3, nf).';
end
z = (z.^(1/degree)) .* reshape(exp(2i * pi .* (0:(degree-1)) / degree), 1, 1, degree);

% stepsize for gradient flow
max_dt = 0.5 * (mean(meshData.edgeLengths) - std(meshData.edgeLengths));

[lb, ub] = bounds(meshData.verts, 1);
curveLength = 0.25 * max(ub - lb);

% number of steps to flow
nsteps = round(curveLength / max_dt);


% initialize paths in interiors of triangles, weighted by area
f_idx = randsample(meshData.nf, ncurves, true, meshData.areas);
bary_local = rand(ncurves, 3);
bary_local = bary_local ./ sum(bary_local, 2);
bary = sparse(faces(f_idx, :), repmat((1:ncurves).', 1, 3), bary_local, nv, ncurves);

vel = randn(ncurves, 1) + 1i * randn(ncurves, 1); % initial direction

paths = nan(ncurves, 3, nsteps);
paths(:, :, 1) = bary' * verts;

good_idx = (1:ncurves).';

% Trace integral curves of directional field
for i=2:nsteps
    % Interpolate field component
    vel = interpVelocity(bary, f_idx, vel);
    
    % Normalize velocity
    vel = vel ./ abs(vel);

    % Calculate step size
    [bary, f_idx, vel, still_good_idx] = maxstep(bary, f_idx, vel);
    good_idx = good_idx(still_good_idx);
    
    % Add 3D points to paths
    paths(good_idx, :, i) = bary' * verts;
end

paths = [permute(paths, [3 1 2]); nan(1, ncurves, 3)];
cmap = cbrewer('Paired', 8);
color_ix = repmat(randi(8, [size(paths, 1), 1]), 1, size(paths, 2)).';
colors = reshape(cmap(color_ix, :), [], 3);
hdl.FaceVertexCData = colors;
hdl = patch(paths(:, :, 1), paths(:, :, 2), paths(:, :, 3), LineWidth=2, ...
    EdgeColor="interp", FaceVertexCData=colors);

% Returns interpolation of the directional field component best matching old velocity
function vel_new = interpVelocity(bary, f_idx, vel)
    nheads = size(f_idx, 1);
    bary_local = full(bary(sub2ind([nv nheads], faces(f_idx, :), repmat((1:nheads).', 1, 3))));
    
    zCandidates = z(f_idx, :, :);
    [~, idx] = max(real(conj(vel) .* zCandidates), [], 3, 'linear');
    vel_new = reshape(zCandidates(idx), nheads, []);
    if size(vel_new, 2) > 1
        vel_new = dot(bary_local, vel_new, 2);
    end
end

function [bary_new, f_idx_new, vel_new, good_idx] = maxstep(bary, f_idx, vel)
    nheads = size(bary, 2);

    bary_local = full(bary(sub2ind([nv nheads], faces(f_idx, :), repmat((1:nheads).', 1, 3))));
    d_bary = G_intr' * sparse(2 * f_idx + (-1:0), repmat(1:nheads, 1, 2), [real(vel) imag(vel)], 2*nf, nheads);
    d_bary = reshape(nonzeros(d_bary), 3, nheads).';

    % Calculate maximum step before any barycentric coordinate becomes zero
    dt = -bary_local ./ d_bary;
    dt(d_bary >= 0) = inf;
    [dt, tight_vtx] = min(dt, [], 2);
    transition_mask = dt < max_dt;
    dt(~transition_mask) = max_dt;

    % Update barycentric coordinates
    bary_local_new = bary_local + dt .* d_bary;
    bary_new = sparse(faces(f_idx, :), repmat((1:nheads).', 1, 3), bary_local_new, nv, nheads);
    
    % hop over edge opposite vertex associated to barycentric coordinate
    % that becomes zero
    tight_edge = sub2ind([nf, 3], f_idx, mod(tight_vtx, 3) + 1);
    f_idx_new = f_idx;
    f_idx_new(transition_mask) = face2face(tight_edge(transition_mask));

    % cut off curves that have gone off a boundary
    good_idx = find(~isnan(f_idx_new));
    f_idx_new = f_idx_new(good_idx);
    bary_new = bary_new(:, good_idx);
    nheads = size(good_idx, 1);

    % parallel-transport current velocity to new face
    % transport coeff is 1 for any face to itself
    transport_coeff = full(reshape(f2fTransport(sub2ind([nf nf], f_idx_new, f_idx(good_idx))), nheads, 1));
    vel_new = transport_coeff .* vel(good_idx);
end

end