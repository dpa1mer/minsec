function faces_experiment_figures(datadir)

for j = 1:5
    for deg = [2 4]
        load(fullfile(datadir, sprintf("face%d_deg%d.mat", j, deg)), "face", "z");
        face = process_surface(face.verts, face.faces, Flat=false);

        render_surface(face, 0, 10, 10, 20); hold on;
        plot_integral_curves(face, z, deg, 1000);
        plot_singularities(face, z, deg, 0.02);
        colormap(cbrewer('RdBu', 5));
        exportgraphics(gca, fullfile(datadir, sprintf("face%d_deg%d_minsec.png", j, deg)));
        close all;
    
        z_mbo = mbo_intrinsic(face, .1, deg);
        render_surface(face, 0, 10, 10, 20); hold on;
        plot_integral_curves(face, z_mbo, deg, 1000);
        plot_singularities(face, z_mbo, deg, 0.02);
        colormap(cbrewer('RdBu', 5));
        exportgraphics(gca, fullfile(datadir, sprintf("face%d_deg%d_mbo_0.1.png", j, deg)));
        close all;
    
        z_mbo = mbo_intrinsic(face, .01, deg);
        render_surface(face, 0, 10, 10, 20); hold on;
        plot_integral_curves(face, z_mbo, deg, 1000);
        plot_singularities(face, z_mbo, deg, 0.02);
        colormap(cbrewer('RdBu', 5));
        exportgraphics(gca, fullfile(datadir, sprintf("face%d_deg%d_mbo_0.01.png", j, deg)));
        close all;

        z_knoppel_ours = mbo_intrinsic(face, -1, deg);
        render_surface(face, 0, 10, 10, 20); hold on;
        plot_integral_curves(face, z_knoppel_ours, deg, 1000);
        plot_singularities(face, z_knoppel_ours, deg, 0.02);
        colormap(cbrewer('RdBu', 5));
        exportgraphics(gca, fullfile(datadir, sprintf("face%d_deg%d_knoppel_ours.png", j, deg)));
        close all;

        % Geometry Central aligns to normals, whereas we align to tangents.
        % So rotate directions by \pi/2.
        z_knoppel_geomcentral = compute_field_geometry_central(face, deg) .* (1i).^deg;
        render_surface(face, 0, 10, 10, 20); hold on;
        plot_integral_curves(face, z_knoppel_geomcentral, deg, 1000);
        plot_singularities(face, z_knoppel_geomcentral, deg, 0.02);
        colormap(cbrewer('RdBu', 5));
        exportgraphics(gca, fullfile(datadir, sprintf("face%d_deg%d_knoppel_geomcentral.png", j, deg)));
        close all;
    end
end

end