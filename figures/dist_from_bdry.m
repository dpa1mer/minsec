function dist = dist_from_bdry(meshData)
% DIST_FROM_BDRY - Compute approximate geodesic distances from the mesh
% boundary using the heat method.
%   dist = DIST_FROM_BDRY(meshData)

assert(meshData.nb > 0);

nv = meshData.nv;
nf = meshData.nf;

bdry_ind = zeros(nv, 1);
bdry_ind(meshData.bdryIdx) = 1;

dt = sqrt(mean(meshData.edgeLengths));

heat = (meshData.star0 + dt * meshData.L) \ (meshData.star0 * bdry_ind);

G = covariant_ops(meshData, 0);
A = spdiags(meshData.areas, 0, nf, nf);

Div = G' * kron(A, speye(2));

V = reshape(G * heat, 2, nf);
V = V ./ vecnorm(V, 2, 1);

L = Div * G;
dist_rhs = Div * V(:);
L_int = L(meshData.intIdx, meshData.intIdx);
dist = zeros(nv, 1);
dist(meshData.intIdx) = -L_int \ dist_rhs(meshData.intIdx);

end