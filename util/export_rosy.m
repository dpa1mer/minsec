function export_rosy(filename, meshData, z, n_components, output_metadata, normalize)

if nargin < 6
    normalize = true;
end

if nargin < 5
    output_metadata = true;
end

if nargin < 4
    n_components = 1;
end
assert(n_components >= 1 && n_components <= 4);

if size(z, 1) == meshData.nv
    [~, AvgV2F] = covariant_ops(meshData, 4);
    z = AvgV2F * z;
end
assert(size(z, 1) == meshData.nf)

if normalize
    z = z ./ vecnorm(z, 2, 2);
end

% Pick one frame component; it does not matter which
z = (z.^(1/4));

z = reshape(z.' .* exp((2 * pi * 1i / 4) .* (0:(n_components-1)).'), 1, n_components * meshData.nf);

% Convert to ambient basis
z = reshape(pagemtimes(meshData.faceBasis, ...
                       reshape([real(z); imag(z)], 2, n_components, meshData.nf)), ...
            3 * n_components, meshData.nf).';

if output_metadata
    f = fopen(filename, "w");
    fprintf(f, "%d\n", meshData.nf);
    fprintf(f, "%d\n", 4);
    fclose(f);
end
writematrix(z, filename, Delimiter=' ', WriteMode='append', FileType='text');

end