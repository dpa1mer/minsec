function [h_lift, h_quiv] = plot_section_annotations(meshData, z, idx, shift, scale, degree, radius)

if nargin < 7
    radius = 300;
end

if nargin < 6
    degree = 1;
end

if nargin < 5
    scale = 2 * pi;
end

if nargin < 4
    shift = 0.5;
end

z = z(idx);
base = repmat(meshData.verts(idx, :), degree, 1);

height = angle(z) ./ (2*pi);
height = height + shift;
height = mod(height, 1);
height = reshape((height + (0:(degree - 1))) * scale, [], 1);
h_lift = scatter3(base(:, 1), base(:, 2), height, radius, 'k.');

v = z.^(1/degree);
v = v .* exp(2i * pi .* (0:(degree-1)) / degree);
v = [real(v(:)), imag(v(:))];
v = mean(meshData.edgeLengths) * v;
h_quiv = quiver(base(:, 1), base(:, 2), v(:, 1), v(:, 2), 'k', 'LineWidth', 2, 'AutoScale', 'off');

end