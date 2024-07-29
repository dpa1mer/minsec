function meshData = process_surface(verts, faces, options)

arguments
    verts (:, 3) double
    faces (:, 3) double {mustBeInteger, mustBeGreaterThan(faces, 0)}

    options.Flat          (1, 1) logical = false
    options.NormalizeArea (1, 1) double
end

% Get rid of hanging triangles and re-index vertices
faces = flip_ears(verts, faces, 'FlipAndClip', true);
[used_vert_ix, ~, faces] = unique(faces);
faces = reshape(faces, [], 3);
verts = verts(used_vert_ix, :);

% Center and rescale so that total area is 1
if isfield(options, 'NormalizeArea')
    M = massmatrix(verts, faces, 'barycentric');
    center = sum(M * verts, 1) / sum(M, 'all');
    verts = verts - center;
    areas = 0.5 * vecnorm(cross(verts(faces(:, 1), :) - verts(faces(:, 2), :), ...
                                verts(faces(:, 3), :) - verts(faces(:, 2), :)), 2, 2);
    totalArea = sum(areas);
    verts = verts .* sqrt(options.NormalizeArea / totalArea);
end

tri = triangulation(faces, verts);
meshData.tri = tri;
meshData.faceNormals = tri.faceNormal;
meshData.vertNormals = tri.vertexNormal;

meshData.verts = verts; meshData.faces = faces;
meshData.faceOrientedEdges = reshape(faces(:, [1 2; 2 3; 3 1]), [], 2);
meshData.faceOrientation = reshape(sign(meshData.faceOrientedEdges(:, 2) - meshData.faceOrientedEdges(:, 1)), [], 3);
faceEdges = sort(meshData.faceOrientedEdges, 2);
[meshData.edges, ~, meshData.face2edge] = unique(faceEdges, 'rows');
meshData.face2edge = reshape(meshData.face2edge, [], 3);
ne = size(meshData.edges, 1); meshData.ne = ne;
nv = size(verts, 1); meshData.nv = nv;
nf = size(faces, 1); meshData.nf = nf;

meshData.edgeLengths = vecnorm(meshData.verts(meshData.edges(:, 2), :) - meshData.verts(meshData.edges(:, 1), :), 2, 2);
meshData.faceCenters = (1/3) * squeeze(sum(reshape(verts(faces, :), nf, 3, 3), 2));

meshData.edgeTangents = meshData.verts(meshData.edges(:, 2), :) - meshData.verts(meshData.edges(:, 1), :);
meshData.edgeTangents = meshData.edgeTangents ./ vecnorm(meshData.edgeTangents, 2, 2);

meshData.edgeMidpoints = (meshData.verts(meshData.edges(:, 1), :) + meshData.verts(meshData.edges(:, 2), :)) / 2;

%% Construct Hodge Stars
vertOppEdges = meshData.face2edge(:, [2 3 1]);
eij = reshape(verts(faces(:, [3 1 2]), :) - verts(faces, :), nf, 3, 3);
eik = reshape(verts(faces(:, [2 3 1]), :) - verts(faces, :), nf, 3, 3);
meshData.faceVertAngles = acos((dot(eij, eik, 3) ./ vecnorm(eij, 2, 3)) ./ vecnorm(eik, 2, 3));
meshData.vertCotans = cot(meshData.faceVertAngles);

Lij = meshData.vertCotans / 2;
L = sparse(meshData.edges(vertOppEdges, 1), meshData.edges(vertOppEdges, 2), Lij, nv, nv);
meshData.star1 = spdiags(nonzeros(L'), 0, ne, ne);
meshData.star1dual = spdiags(diag(meshData.star1).^(-1), 0, ne, ne);

meshData.areas = 0.5 * vecnorm(cross(verts(faces(:, 1), :) - verts(faces(:, 2), :), verts(faces(:, 3), :) - verts(faces(:, 2), :)), 2, 2);
meshData.star2 = spdiags(meshData.areas.^(-1), 0, nf, nf);

meshData.star0 = massmatrix(verts, faces, 'full');
meshData.star0lump = massmatrix(verts, faces, 'barycentric');

%% Construct Laplacians
meshData.d0 = sparse(repmat((1:ne).', 1, 2), meshData.edges, repmat([-1 1], ne, 1), ne, nv);
meshData.d1 = sparse(repmat((1:nf).', 1, 3), meshData.face2edge, meshData.faceOrientation, nf, ne);
meshData.L = meshData.d0.' * meshData.star1 * meshData.d0;
meshData.Ldual = meshData.d1 * meshData.star1dual * meshData.d1.';

%% Construct boundary matrices
meshData.bdryEdges = freeBoundary(tri);
nb = size(meshData.bdryEdges, 1); meshData.nb = nb;
meshData.ni = nv - nb;

if nb > 0
    meshData.bdryIdx = meshData.bdryEdges(:, 1);
    meshData.intIdx = setdiff(1:meshData.nv, meshData.bdryIdx);

    [bdryForward, bdryEdgeIdxF] = ismember(meshData.bdryEdges, meshData.edges, 'rows');
    [bdryBackward, bdryEdgeIdxB] = ismember(fliplr(meshData.bdryEdges), meshData.edges, 'rows');
    meshData.bdryEdgeIdx = bdryEdgeIdxF + bdryEdgeIdxB;
    meshData.intEdgeIdx = setdiff((1:meshData.ne).', meshData.bdryEdgeIdx);
else
    meshData.bdryIdx = [];
    meshData.intIdx = (1:meshData.nv).';
    meshData.bdryEdgeIdx = [];
    meshData.intEdgeIdx = (1:meshData.ne).';
end

%% Face-edge Levi-Civita connection
meshData.faceEdgeVectors = permute(eik, [3 2 1]); % 3(D) x 3 x nf
meshData.faceEdgeTangents = meshData.faceEdgeVectors ./ vecnorm(meshData.faceEdgeVectors, 2, 1);
meshData.faceEdgeNormals = cross(repmat(reshape(meshData.faceNormals.', 3, 1, nf), 1, 3, 1), ...
                                 meshData.faceEdgeTangents, 1);
meshData.faceEdgeBasis = cat(2, reshape(meshData.faceEdgeTangents, 3, 1, 3, nf), ...
                                reshape(meshData.faceEdgeNormals, 3, 1, 3, nf));
meshData.faceBasis = squeeze(meshData.faceEdgeBasis(:, :, 1, :));
if options.Flat
    meshData.faceBasis = repmat([1 0; 0 1; 0 0], 1, 1, nf);
end

relAngles = [zeros(nf, 1), -pi + meshData.faceVertAngles(:, 2), - pi - meshData.faceVertAngles(:, 1)];
meshData.f2eTransportAngle = relAngles + angle(meshData.faceOrientation);
meshData.f2eTransport = meshData.faceOrientation .* exp(1i .* relAngles);
if options.Flat
    meshData.f2eTransportAngle = zeros(nf, 3);
    meshData.f2eTransport = ones(nf, 3);
end

meshData.angleDefect = (2*pi) - accumarray(faces(:), meshData.faceVertAngles(:), [nv 1]);
if nb > 0
    meshData.angleDefect(meshData.bdryIdx) = meshData.angleDefect(meshData.bdryIdx) - pi;
end

%% Face-vertex connection

% For vertex frames, pick the edge of lowest index, counting ingoing edges
% as negative (this choice is arbitrary).
meshData.vertBasisEdge = accumarray(meshData.edges(:), reshape((1:ne).' .* [1, -1], [], 1), [nv 1], @min);
vertBasisX = meshData.edgeTangents(abs(meshData.vertBasisEdge), :) .* ...
             sign(meshData.vertBasisEdge);
vertBasisX = vertBasisX - dot(vertBasisX, meshData.vertNormals, 2) .* meshData.vertNormals;
vertBasisX = vertBasisX ./ vecnorm(vertBasisX, 2, 2);
if options.Flat
    vertBasisX = repmat([1 0 0], nv, 1);
end
vertBasisY = cross(meshData.vertNormals, vertBasisX, 2);
meshData.vertBasis = permute(cat(3, vertBasisX, vertBasisY), [2 3 1]); % 3 x 2 x nv

[~, ~, f2vRotationExtrinsic] = principal_rotation( ...
    repmat(reshape(meshData.faceNormals.', 3, 1, nf), 1, 3, 1), ...
    permute(reshape(meshData.vertNormals(meshData.faces, :), nf, 3, 3), [3 2 1]));
v2fRotation = pagemtimes( ...
    pagemtimes(f2vRotationExtrinsic, reshape(meshData.faceBasis, 3, 2, 1, nf)), 'transpose', ...
    reshape(meshData.vertBasis(:, :, meshData.faces.'), 3, 2, 3, nf), 'none');
meshData.v2fTransport = squeeze(v2fRotation(1, 1, :, :) + 1i .* v2fRotation(2, 1, :, :));
meshData.v2fTransportAngle = angle(meshData.v2fTransport);

%% Vertex-vertex connection
[~, ~, v2vRotationExtrinsic] = principal_rotation( ...
    meshData.vertNormals(meshData.edges(:, 1), :).', ...
    meshData.vertNormals(meshData.edges(:, 2), :).');
v2vRotation = pagemtimes( ...
    meshData.vertBasis(:, :, meshData.edges(:, 2)), 'transpose', ...
    pagemtimes(v2vRotationExtrinsic, meshData.vertBasis(:, :, meshData.edges(:, 1))), 'none');
meshData.v2vTransport = reshape(v2vRotation(1, 1, :) + 1i .* v2vRotation(2, 1, :), [], 1);
meshData.v2vTransportAngle = angle(meshData.v2vTransport);

end