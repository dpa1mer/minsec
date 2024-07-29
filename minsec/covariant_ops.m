function [V2F, AvgV2F, AvgF2V, F2C] = covariant_ops(meshData, k)

K = length(k);
nf = meshData.nf;
nv = meshData.nv;

% Parallel transport values at vertices to corners
findex = reshape((1:(3 * nf * K)), 3, nf, K);
vindex = reshape((0:K-1) * nv, 1, 1, K) + meshData.faces.';
v2f = exp(1i .* reshape(k, 1, 1, K) .* meshData.v2fTransportAngle);
V2F = sparse(findex(:), vindex(:), v2f(:), 3 * nf * K, nv * K);

% Explode face values to corners
F2C = kron(speye(nf), ones(3, 1));

% Averaging operator - maps field values at verts to faces covariantly
sumPerFace = kron(speye(K), F2C');
AvgV2F = (1/3) * sumPerFace * V2F;
AvgF2V = (1/3) * kron(speye(K), meshData.star0lump^(-1)) * V2F' * sumPerFace' * kron(speye(K), meshData.star2^(-1));


end