function field = compute_field_directional(meshData, degree, bdryVals)

persistent mh
if ~(isa(mh,'matlab.mex.MexHost') && isvalid(mh))
    mh = mexhost;
end

% The mex function returns field representatives in extrinsic coordinates
field = feval(mh, "compute_field_directional_mex", meshData, degree, bdryVals);

% Convert to intrinsic power field
field = squeeze(pagemtimes(meshData.vertBasis, 'transpose', ...
                           reshape(field(:, 1:3).', 3, 1, meshData.nv), 'none')).';
field = (field(:, 1) + 1i * field(:, 2)).^degree;

if any(isnan(field))
    warning("Field contains NaN values!")
end

end