function qvr = fancy_quiver(points, vectors, colorMap, cutoff, varargin)

if nargin < 3
    colorMap = inferno;
end

norms = vecnorm(vectors, 2, 2);
if nargin < 4
    cutoff = mean(norms);% + std(norms);
end

ix = find(norms >= cutoff);
if numel(ix) <= 1
    return;
end

qvr = quiver3(points(ix,1), points(ix, 2), points(ix, 3), vectors(ix, 1), vectors(ix, 2), vectors(ix, 3), 'LineWidth', 1.5, 'Color', 'k', varargin{:});

% Color a la
% https://stackoverflow.com/questions/29632430/quiver3-arrow-color-corresponding-to-magnitude
colorField = norms(ix, :);
% [~, perm] = sort(colorField);
% colorField(perm) = 1:length(perm);
[~, ~, colorIdx] = histcounts(colorField, size(colorMap, 1));
cmap = uint8(ind2rgb(colorIdx(:), colorMap) * 255);
cmap(:, :, 4) = 255; % Opaque
cmap = permute(repmat(cmap, [1 2 1]), [2 1 3]);
qvr.Tail.ColorData = reshape(cmap, [], 4)';
qvr.Tail.ColorBinding = 'interpolated';

end