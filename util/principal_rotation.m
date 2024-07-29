function [axis, angle, rotm] = principal_rotation(vec1, vec2)

assert(size(vec1, 1) == 3 && size(vec2, 1) == 3);
assert(all(size(vec1) == size(vec2)));
assert(all(abs(vecnorm(vec1, 2, 1) - 1) < sqrt(eps), 'all'));
assert(all(abs(vecnorm(vec2, 2, 1) - 1) < sqrt(eps), 'all'));

orig_size = size(vec1);

vec1 = reshape(vec1, 3, []); vec2 = reshape(vec2, 3, []);

axis = cross(vec1, vec2, 1);
angle = atan2(vecnorm(axis, 2, 1), dot(vec1, vec2, 1));
axis = axis ./ vecnorm(axis, 2, 1);

% Fix NaNs when angle is zero
axis(isnan(axis)) = 0;

axis = reshape(axis, orig_size);
angle = reshape(angle, [1 orig_size(2:end)]);

if nargout > 2
    rotm = rotate_by( ...
        repmat(reshape(axis, [3 1 orig_size(2:end)]), 1, 3), ...
        repmat(reshape(angle, [1 1 orig_size(2:end)]), 1, 3), ...
        repmat(eye(3), [1 1 orig_size(2:end)]));
end

end