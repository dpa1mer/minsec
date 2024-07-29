function plot_singular_measure(meshData, C)

ne = meshData.ne;
nf = meshData.nf;
faces = reshape(1:(3*nf), 3, nf).';
verts = meshData.verts(meshData.faces.', :);

edgeValues = zeros(ne, 1);
edgeValues(meshData.intEdgeIdx) = C;
vertValues = [1 -1 1; 1 1 -1; -1 1 1] * edgeValues(meshData.face2edge.');

patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', vertValues(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap viridis;

end