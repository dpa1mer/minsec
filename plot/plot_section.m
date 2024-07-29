function ptc = plot_section(meshData, S, ptc, options)

arguments
    meshData struct
    S        double
    ptc = []

    options.BaseSurf      (1, 1) logical = true
    options.Degree        (1, 1) double {mustBeInteger,mustBePositive} = 1
    options.Normalized    (1, 1) logical = true
    options.SubsampleRate (1, 1) double { ...
        mustBePositive, ...
        mustBeLessThanOrEqual(options.SubsampleRate, 1)} = 1
    options.NormalShift (1, 1) double = 0.1
end

if ndims(S) == 4
    S = squeeze(mean(S, 2));
end

if options.Normalized
    S = S ./ max(vecnorm(S, 2, 1), [], 3);
end

persistent cmap cmap_size;
if isempty(cmap)
    cmap_size = 1024;
    cmap = cbrewer('YlGnBu', cmap_size);
end

% Repeat values degree times
S = repmat(S, 1, 1, options.Degree);

% S is a 1-form/vector field on the circle bundle over the surface,
% represented per face, i.e., of size 3 x nf x K
nf = meshData.nf;
assert(size(S, 2) == nf);
K = size(S, 3);

% subsample faces
sampledFaceIdx = find(rand(nf, 1) <= options.SubsampleRate);%randi(nf, [nSamples, 1]);
nSamples = numel(sampledFaceIdx);
sampledFaceCenters = meshData.faceCenters(sampledFaceIdx, :);

% offset a little in the normal direction to prevent clipping
sampledFaceCenters = sampledFaceCenters + ...
    options.NormalShift * mean(meshData.edgeLengths) * meshData.faceNormals(sampledFaceIdx, :);

faceBasis = meshData.faceBasis(:, :, sampledFaceIdx);
S = S(:, sampledFaceIdx, :);

normS = vecnorm(S, 2, 1);
plotCoordsComplex = permute(normS .* exp(2i * pi .* reshape(0:(K-1), 1, 1, K) / K), [1 3 2]);
plotCoords2D = [real(plotCoordsComplex); imag(plotCoordsComplex)];

scale = (0.3 / options.SubsampleRate) * mean(meshData.edgeLengths);
plotCoords3D = scale * pagemtimes(faceBasis, plotCoords2D) + reshape(sampledFaceCenters.', 3, 1, nSamples);

% Add vertex at origin for each plot
plotVerts = [reshape(sampledFaceCenters.', 3, 1, nSamples), plotCoords3D];
plotVerts = reshape(plotVerts, 3, (K + 1) * nSamples).';

plotFaces = [(2:(K+1)); circshift(2:(K+1), -1, 2); ones(1, K)] + ...
            (K + 1) * reshape((0:(nSamples - 1)), 1, 1, nSamples);
plotFaces = reshape(plotFaces, 3, K * nSamples).';

colors = reshape([zeros(1, nSamples); squeeze(normS).'], (K + 1) * nSamples, 1);
cmin = min(colors);
cmax = max(colors);
cidx = min(cmap_size, round((cmap_size - 1) * (colors - cmin) / (cmax - cmin)) + 1);
colors = squeeze(ind2rgb(cidx, cmap));

if isempty(ptc)
    if options.BaseSurf
        trisurf(meshData.tri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on;
    end
    ptc = patch('Faces', plotFaces, 'Vertices', plotVerts, 'EdgeColor', 'none', ...
                'FaceColor', 'flat', 'FaceVertexCData', colors);
else % Update existing plot
    ptc.Faces = plotFaces;
    ptc.Vertices = plotVerts;
    ptc.FaceVertexCData = colors;
end

end