function qvr = plot_intrinsic_field(meshData, v, degree, options)

arguments
    meshData struct
    v        double
    degree  (1, 1) double {mustBePositive,mustBeInteger} = 1

    options.ColorMap          (:, :) double = inferno
    options.BackgroundSurface (1, 1) logical = true
    options.NormalShift       (1, 1) double = 0
end

if size(v, 1) == meshData.nf
    pts = meshData.faceCenters + options.NormalShift * mean(meshData.edgeLengths) *  meshData.faceNormals;
    frames = meshData.faceBasis;
elseif size(v, 1) == meshData.nv
    pts = meshData.verts + options.NormalShift * mean(meshData.edgeLengths) *  meshData.vertNormals;
    frames = meshData.vertBasis;
else
    error('Vector field is wrong size');
end

nPts = size(pts, 1);
pts = repmat(pts, degree, 1);
frames = repmat(frames, 1, 1, degree);

assert(size(v, 2) == 1);
v = v.^(1/degree);
v = v .* exp(2i * pi .* (0:(degree-1)) / degree);
v = [real(v(:)).'; imag(v(:)).'];

v = squeeze(pagemtimes(frames, reshape(v, 2, 1, degree * nPts))).';
if options.BackgroundSurface
    trisurf(meshData.tri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on;
end
qvr = fancy_quiver(pts, v, options.ColorMap, 0, ShowArrowHead=false);

end