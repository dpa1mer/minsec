function [bdryVertAngles, bdryVertPhases, BdryEdgeIntOp, bdryEdgeIntVal] = field_bdry_conds(meshData, degree)

nf = meshData.nf;
nb = meshData.nb;

[isBdryEdge, bdryF2E] = ismember(meshData.face2edge, meshData.bdryEdgeIdx);

edgeProj = pagemtimes(meshData.faceBasis, 'transpose', meshData.faceEdgeVectors, 'none');
edgeProj = reshape(permute(edgeProj, [2 3 1]), 3 * nf, 2);
bdryProj = edgeProj(isBdryEdge.', :);

intrinsicEdgeTangents = reshape(pagemtimes(meshData.faceBasis, 'transpose', ...
                                           meshData.faceEdgeTangents, 'none'), 2, 3 * nf);
bdryFaceTangents = intrinsicEdgeTangents(:, isBdryEdge.');

% This construction assumes that the boundary and faces are oriented
% consistently (counter-clockwise)
bdryFaces = find(any(isBdryEdge, 2));
bdryF2E = max(bdryF2E(bdryFaces, :), [], 2);
BdryEdgeIntOp = sparse(repmat(bdryF2E, 1, 2), 2 * (bdryFaces - 1) + (1:2), bdryProj, nb, nf * 2);

[V2F, ~, AvgF2V] = covariant_ops(meshData, degree);
bdryFacePhases = ((bdryFaceTangents(1, :) + 1i * bdryFaceTangents(2, :)).^degree).';
bdryVertPhases = AvgF2V(meshData.bdryIdx, bdryFaces) * bdryFacePhases;
bdryVertPhases = bdryVertPhases ./ abs(bdryVertPhases);
bdryVertAngles = angle(bdryVertPhases);

bdryAngleDiff = kron(speye(nf), [-1 1 0; 0 -1 1; 1 0 -1]) * ...
                angle(V2F(:, meshData.bdryIdx) * bdryVertPhases);
bdryAngleDiff = bdryAngleDiff(isBdryEdge.');

bdryEdgeIntVal = zeros(nb, 1);
bdryEdgeIntVal(bdryF2E) = -bdryAngleDiff / (2 * pi);
bdryEdgeIntVal = bdryEdgeIntVal - round(bdryEdgeIntVal);

end