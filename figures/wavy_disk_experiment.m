function wavy_disk_experiment(outdir, filename, eccentricity)

[raw_verts, faces] = readOBJ("models/disk-mini.obj");
[raw_verts, faces] = loop(raw_verts, faces);
[raw_verts, faces] = loop(raw_verts, faces);

for t = 0:0.4:0.8
    verts(:, 1) = raw_verts(:, 1);
    verts(:, 2) = raw_verts(:, 2) * eccentricity;
    verts(:, 3) = t * (.1 * sin(4*pi*raw_verts(:, 1)) - raw_verts(:, 2).^2);
    wavy_disk = process_surface(verts, faces, Flat=false, NormalizeArea=1);
    % figure; trisurf(wavy_disk.tri); view(3); axis image vis3d off;
    
    data_file = fullfile(outdir, sprintf("%s_%.1f.mat", filename, t));
    if isfile(data_file)
        load(data_file, "wavy_disk", "S", "C", "f", "z", "admmData");
        wavy_disk = process_surface(wavy_disk.verts, wavy_disk.faces, Flat=false);
    else
        [S, C, f, z, admmData] = minsec(wavy_disk, 64, .1, 4, .1);
        save(data_file, "wavy_disk", "S", "C", "f", "z", "admmData");
    end
    
    render_surface(wavy_disk, 30, 30, 10, 20);
    hold on; plot_integral_curves(wavy_disk, z, 4, 2000);
    plot_singularities(wavy_disk, z, 4, 0.015);
    colormap(cbrewer('RdBu', 5));
    exportgraphics(gca, fullfile(outdir, sprintf("%s_%.1f.png", filename, t)));
end

end