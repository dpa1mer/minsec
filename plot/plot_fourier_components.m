function plot_fourier_components(meshData, f, K, options)

arguments
    meshData struct
    f           (:, :) double
    K           (1, 1) double {mustBeInteger, mustBePositive} = size(f, 2)

    options.DisplayType (1, 1) string {mustBeMember(options.DisplayType, ["complex", "norm", "angle"])} = "complex"
    options.kthRoot     (1, 1) logical = false
    options.Flat        (1, 1) logical = false
    options.Separate    (1, 1) logical = false
end

nv = meshData.nv;
assert(size(f, 1) == nv);
assert(K <= size(f, 2));

m = floor(sqrt(K));
n = ceil(K / m);
if ~options.Separate
    fig = figure(Color='w'); tl = tiledlayout(m, n, 'TileSpacing', 'compact', 'Padding', 'compact');
end
for k = 1:K
    if options.Separate
        figure(Color='w');
    else
        tiles(k) = nexttile(k);
    end

    if k == 1
        trisurf(meshData.tri, 'FaceVertexCData', f(:, 1), 'FaceColor', 'interp', 'EdgeColor', 'none');
        colormap gray;
    else
        vals = f(:, k);
        if options.kthRoot
            vals = vals .^ (1/(k - 1));
        end

        if options.DisplayType == "angle"
            norms = ones(nv, 1);
        else
            norms = abs(vals) / max(abs(vals));
        end

        angs = mod(angle(vals), 2 * pi) / (2 * pi);

        if options.DisplayType == "norm"
            trisurf(meshData.tri, 'FaceVertexCData', norms, 'FaceColor', 'interp', 'EdgeColor', 'none');
            colormap viridis;
        else
            colors = hsv2rgb([angs, ones(nv, 1), norms]);
            trisurf(meshData.tri, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'EdgeColor', 'none');
        end
    end
    if options.Flat
        view(2); axis image off;
    else
        view(3); axis image vis3d off;
    end
end

if ~options.Separate
    Link = linkprop(tiles,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(fig, 'StoreTheLink', Link);
end

end