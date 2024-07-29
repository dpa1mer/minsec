function [V, D, F] = quad_mesh_directional(meshData, z, options)
% QUAD_MESH_DIRECTIONAL Generate a quad mesh from a cross field using the Directional
% library.
%   [V, D, F] = QUAD_MESH(meshData, z, options)
%
%   V: vertices of the generated mesh
%   D: face degrees
%   F: faces as indices into V


arguments
    meshData    (1, 1) struct
    z           (:, 1) double
    
    options.normalize       (1, 1) logical  = true
    options.verbose         (1, 1) logical  = false
    options.seamless        (1, 1) logical  = true
    options.round_seams     (1, 1) logical  = false
    options.length_ratio    (1, 1) double   = 0.02
end

persistent mh
if ~(isa(mh,'matlab.mex.MexHost') && isvalid(mh))
    mh = mexhost;
end

nv = meshData.nv;
nf = meshData.nf;

if options.normalize
    z = z ./ abs(z);
end

if size(z, 1) == nv
    [~, AvgV2F] = covariant_ops(meshData, 4);
    z = AvgV2F * z;
end

z = reshape((z.').^(1/4) .* [1; 1i; -1; -1i], 1, 4, nf);
field = pagemtimes(meshData.faceBasis, [real(z); imag(z)]);
field = reshape(field, 3 * 4, nf).';

[V, D, F] = feval(mh, "quad_mesh_directional_mex", meshData, field, options);

F = double(F);
F(1:size(F, 2) > D) = nan;
D = double(D);

end