function boundary_opening_experiment(outdir)

assert(isfolder(outdir));

[verts, faces] = readOBJ('../models/sphere.obj');
for offset = [-0.9, -0.5, 0]
    [verts_cut, faces_cut] = half_space_intersect(verts, faces, [0 0 offset], [0 0 1], Cap=false);
    sphere_cut = process_surface(verts_cut, faces_cut);

    [S, C, f, z, admmData] = minsec(sphere_cut, 64, .1, 1, .1);
    save(fullfile(outdir, sprintf("sphere_cut_offset%0.2f_deg1.mat", offset)), "sphere_cut", "S", "C", "f", "z", "admmData");

    render_surface(sphere_cut, 0, 10, 0, 20);
    hold on; plot_section(sphere_cut, S, BaseSurf=false);
    [~, AvgV2F] = covariant_ops(sphere_cut, 1);
    qvr = plot_intrinsic_field(sphere_cut, AvgV2F * z, 1, NormalShift=0.1, ColorMap=[1 0 0]);
    material dull;
    exportgraphics(gca, fullfile(outdir, sprintf("sphere_cut_offset%0.2f_deg1_sec.png", offset)), Resolution=300);

    render_surface(sphere_cut, 0, 10, 0, 20);
    hold on; plot_integral_curves(sphere_cut, z, 1, 2000);
    material dull;
    exportgraphics(gca, fullfile(outdir, sprintf("sphere_cut_offset%0.2f_deg1_int.png", offset)), Resolution=300);

    close all;
end

end