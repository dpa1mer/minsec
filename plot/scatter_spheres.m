function scatter_spheres(centers, colors, radii, n_verts)

if nargin < 4
    n_verts = 20;
end

persistent sph_xyz n_sphere_verts
if isempty(sph_xyz) || isempty(n_sphere_verts) || n_sphere_verts ~= n_verts
    n_sphere_verts = n_verts;
    [sph_x, sph_y, sph_z] = sphere(n_sphere_verts);
    sph_xyz = cat(3, sph_x, sph_y, sph_z);
    sph_xyz = [sph_xyz nan(n_sphere_verts + 1, 1, 3)];
    sph_xyz = reshape(sph_xyz, n_sphere_verts + 1, n_sphere_verts + 2, 1, 3);
end

N = size(centers, 1);
assert(size(centers, 2) == 3)
centers = reshape(centers, 1, 1, N, 3);

if ~isscalar(radii)
    assert(all(size(radii) == [N 1]));
    radii = reshape(radii, 1, 1, N);
end

scatter_xyz = reshape((sph_xyz .* radii) + centers, n_verts + 1, (n_verts + 2) * N, 3);

if size(colors, 1) == 1
    colors = repmat(colors, N, 1);
end

colors_size = size(colors);
assert(N == colors_size(1));
scatter_color = repmat(reshape(colors, [1 1 N colors_size(2:end)]), n_verts + 1, n_verts + 2);
scatter_color = reshape(scatter_color, [n_verts + 1, (n_verts + 2) * N, colors_size(2:end)]);

surf(scatter_xyz(:, :, 1), scatter_xyz(:, :, 2), scatter_xyz(:, :, 3), ...
     scatter_color, 'EdgeColor', 'none', 'FaceColor', 'interp');

end