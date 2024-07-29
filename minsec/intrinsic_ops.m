function [G_intr, M_intr, G_nc, M_nc] = intrinsic_ops(meshData)

nf = meshData.nf;
nv = meshData.nv;
ne = meshData.ne;

% Intrinsic gradient operator
faceEdgeAltitudes = 2 * meshData.areas ./ meshData.edgeLengths(meshData.face2edge);
hatGradients = meshData.faceEdgeNormals ./ reshape(faceEdgeAltitudes.', 1, 3, nf);
hatGradientsIntrinsic = pagemtimes(meshData.faceBasis, 'transpose', hatGradients, 'none');
findex = repmat(reshape(1:2*nf, 2, 1, nf), 1, 3, 1);
fvindex = repmat(reshape([3; 1; 2] + 3 * (0:(nf-1)), 1, 3, nf), 2, 1, 1); % Verts opposite edges
G_intr = sparse(findex(:), fvindex(:), hatGradientsIntrinsic(:), 2 * nf, 3 * nf);

M_intr = kron(spdiags(meshData.areas, 0, nf, nf), [1/6 1/12 1/12; 1/12 1/6 1/12; 1/12 1/12 1/6]);

% Nonconforming intrinsic gradient operator
eindex = repmat(reshape(meshData.face2edge.', 1, 3, nf), 2, 1, 1);
G_nc = sparse(findex(:), eindex(:), -2 * hatGradientsIntrinsic(:), 2 * nf, ne);

edgeAreas = repmat(meshData.areas, 1, 3);
edgeAreas = accumarray(meshData.face2edge(:), edgeAreas(:), [ne 1]) / 3;
M_nc = spdiags(edgeAreas, 0, ne, ne);

end