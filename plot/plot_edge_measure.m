function plot_edge_measure(meshData, C, M)

ne = meshData.ne;
assert(size(C, 1) == ne);

% Sample edge centers in proportion to their weight
sgn = sign(C);

Cabs = abs(C);

ptIdx = randsample(ne, 1000, true, M * Cabs);
pts = meshData.edgeMidpoints(ptIdx, :);

figure; scatter3(pts(:, 1), pts(:, 2), pts(:, 3), 100, sgn(ptIdx), '.');
colormap([1 0 0; 0 0 1]);

end