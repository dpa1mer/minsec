function plt = plot_section_volume(meshData, intensity, plt, options)

arguments
    meshData struct
    intensity double
    plt = []

    options.Scale  (1, 1) double {mustBePositive}               = 1.0
    options.Repeat (1, 1) double {mustBeInteger,mustBePositive} = 1

    % Shift in units of fiber height
    options.Shift  (1, 1) double                                = 0.0

    options.Grid   (:, 1) double {mustBeInteger,mustBePositive} = []
end

repeat = options.Repeat;
fiberScale = options.Scale * repeat * 2 * pi;

nv = meshData.nv;
nf = meshData.nf;

if size(intensity, 1) == 2 || size(intensity, 1) == 3
    intensity = squeeze(vecnorm(intensity, 2, 1));
end

if ~isempty(options.Grid)
    M = options.Grid;
    N = size(intensity, ndims(intensity));
    [grx, gry] = ndgrid(linspace(-1, 1, M));
    grid = cat(3, grx, gry);
    grid = reshape(grid, [], 2);
    tri = triangulation(meshData.faces, meshData.verts(:, 1:2));
    [ix, bary] = pointLocation(tri, grid);
    bary = bary(~isnan(ix), :);
    intensity_grid = zeros(size(grid, 1), 64);
    if ndims(intensity) == 3
        intensity_grid(~isnan(ix), :) = squeeze(sum(bary.' .* intensity(:, ix(~isnan(ix)), :)));
    else
        corner_ix = meshData.faces(ix(~isnan(ix)), :);
        intensity_grid(~isnan(ix), :) = squeeze(sum(bary .* reshape(intensity(corner_ix, :), [size(corner_ix), size(intensity, 2)]), 2));
    end
    intensity_grid = reshape(intensity_grid, M, M, N);
    intensity_grid = repmat(intensity_grid, 1, 1, options.Repeat);
    plt = intensity_grid;
    return;
end

if ndims(intensity) == 3
    assert(size(intensity, 1) == 3);
    assert(size(intensity, 2) == nf);
    intensity = reshape(intensity, 3 * nf, []);
elseif ndims(intensity) == 2
    if size(intensity, 1) == nf
        intensity = repelem(intensity, [3 1]);
    elseif size(intensity, 1) == nv
        intensity = intensity(meshData.faces.', :);
    end
end

N = size(intensity, 2) * repeat;

bundleVerts = reshape(meshData.verts(meshData.faces.', :).' + ...
                      reshape([zeros(2, N); fiberScale * (0:(N-1)) / N], 3, 1, N), 3, []).';
bundleVerts = bundleVerts + [0 0 (options.Shift * 2 * pi * options.Scale)];
bundleFaces = reshape(reshape(1:3*nf, 3, nf) + reshape(3 * nf .* (0:(N-1)), 1, 1, N), 3, []).';

intensity = repmat(intensity(:), repeat, 1);
color = intensity;
alpha = intensity / max(intensity, [], "all");

if isempty(plt)
    plt = patch('Faces', bundleFaces, 'Vertices', bundleVerts, 'FaceVertexCData', color, ...
                'FaceColor', 'interp', 'EdgeColor', 'none', ...
                'FaceAlpha', 'interp', 'AlphaData', 'scaled', 'FaceVertexAlphaData', alpha);
else % Update existing plot
    plt.FaceVertexCData = color;
    plt.FaceVertexAlphaData = alpha;
end

end