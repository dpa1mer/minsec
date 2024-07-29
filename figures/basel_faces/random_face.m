function face = random_face(display, exaggerate)

if nargin < 1
    display = false;
end

if nargin < 2
    exaggerate = 1;
end

persistent faces vertsMean vertsPCA vertsPCAVar
if isempty(faces)
    faces = h5read('model2019_face12.h5', '/shape/representer/cells');
    faces = faces + 1;
    vertsMean = h5read('model2019_face12.h5', '/shape/model/mean');
    vertsPCA = h5read('model2019_face12.h5', '/shape/model/pcaBasis');
    vertsPCAVar = h5read('model2019_face12.h5', '/shape/model/pcaVariance');
end

rot = [1 0 0; 0 0 1; 0 -1 0];

verts = reshape(vertsMean + ((exaggerate * sqrt(vertsPCAVar) .* randn(size(vertsPCAVar)))' * vertsPCA).', 3, []).';
verts = verts * rot;
face = process_surface(verts, faces, NormalizeArea=1);

if display
    plot_face(face);
end

end