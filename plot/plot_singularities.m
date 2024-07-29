function plot_singularities(meshData, z, degree, radius)

if nargin < 4
    radius = 0.02;
end

if(any(isnan(z)))
    error("z has NaN values!");
end

% Normalize
connection = angle(z(meshData.edges(:, 2)) ./ (meshData.v2vTransport.^degree .* z(meshData.edges(:, 1))));

holonomy = round((meshData.d1 * connection) / (2 * pi));
singular_faces = find(holonomy ~= 0);
singular_points = meshData.faceCenters(singular_faces, :);
singular_holos = holonomy(singular_faces);
scatter_spheres(singular_points, singular_holos, radius);

end