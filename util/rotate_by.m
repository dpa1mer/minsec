function vecs = rotate_by(axis, angle, vecs)

assert((size(axis, 1) == 3) && size(vecs, 1) == 3);
assert((size(axis, 3) == size(angle, 3)) && (size(axis, 3) == size(vecs, 3)));

% Rotate vecs by angle around axis using Rodrigues' formula
vecs = cos(angle) .* vecs + ...
       sin(angle) .* cross(axis, vecs, 1) + ...
       (1 - cos(angle)) .* dot(axis, vecs, 1) .* axis;

end