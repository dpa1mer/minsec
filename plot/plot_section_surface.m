function ptc = plot_section_surface(meshData, z, shift, scale, repeat)

if nargin < 5
    repeat = 1;
end

if nargin < 4
    scale = 2 * pi;
end

if nargin < 3
    shift = 0.5;
end

height = angle(z(meshData.faces)).' ./ (2*pi);
height = height + shift;
height = mod(height, 1);

% Deal with triangles that straddle the periodic boundary
idx = any(height > 0.8, 1) & any(height < 0.2, 1);
badheight = height(:, idx);
badheight(badheight < 0.2) = badheight(badheight < 0.2) + 1;
height(:, idx) = badheight;

plotVerts = [meshData.verts(meshData.faces.', 1:2) height(:)];
plotVerts = reshape(plotVerts, [], 1, 3) + reshape([zeros(repeat, 2) (1:repeat).'], 1, repeat, 3);
plotVerts = reshape(plotVerts, [], 3);

heightNormalized = plotVerts(:, 3) - 1;
heightNormalized = 2 * heightNormalized ./ repeat - 1;
plotAlpha = exp(.1 - .1 ./ max(0, 1 - heightNormalized.^2));
assert(all(isfinite(plotAlpha)));

plotVerts(:, 3) = (plotVerts(:, 3) - 1) * scale;

plotFaces = reshape(1:length(plotVerts), 3, []).';

ptc = patch('Faces', plotFaces, 'Vertices', plotVerts, 'EdgeColor', 'none', ...
            'FaceColor', 'interp', 'FaceVertexCData', plotVerts(:, 3), ...
            'FaceAlpha', 'interp', 'FaceVertexAlphaData', plotAlpha, 'AlphaData', 'scaled');
colormap hsv;

end