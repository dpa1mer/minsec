function dataset_figures(datadir)

data_files = dir(fullfile(datadir, "*.mat"));

for j = 1:length(data_files)
    close all;
    the_file = data_files(j);
    filename = fullfile(the_file.folder, the_file.name);

    [~, pure_name, ~] = fileparts(the_file.name);
    minsec_image_name = fullfile(the_file.folder, strcat(pure_name, "_minsec.png"));
    mbo_image_name = fullfile(the_file.folder, strcat(pure_name, "_mbo.png"));
    knoppel_ours_image_name = fullfile(the_file.folder, strcat(pure_name, "_knoppel_ours.png"));
    knoppel_geomcentral_image_name = fullfile(the_file.folder, strcat(pure_name, "_knoppel_geomcentral.png"));

    load(filename, "meshData", "z", "z_mbo", "z_knoppel");
    z = gather(z);

    render_surface(meshData, 90, 10, 10, 20); hold on;
    plot_integral_curves(meshData, z, 4, 1000);
    plot_singularities(meshData, z, 4, 0.015);
    colormap(cbrewer('RdBu', 5));
    exportgraphics(gca, minsec_image_name, Resolution=300);

    render_surface(meshData, 90, 10, 10, 20); hold on;
    plot_integral_curves(meshData, z_mbo, 4, 1000);
    plot_singularities(meshData, z_mbo, 4, 0.015);
    colormap(cbrewer('RdBu', 5));
    exportgraphics(gca, mbo_image_name, Resolution=300);

    render_surface(meshData, 90, 10, 10, 20); hold on;
    plot_integral_curves(meshData, z_knoppel, 4, 1000);
    plot_singularities(meshData, z_knoppel, 4, 0.015);
    colormap(cbrewer('RdBu', 5));
    exportgraphics(gca, knoppel_ours_image_name, Resolution=300);

    z_knoppel_geomcentral = compute_field_geometry_central(meshData, 4);
    render_surface(meshData, 90, 10, 10, 20); hold on;
    plot_integral_curves(meshData, z_knoppel_geomcentral, 4, 1000);
    plot_singularities(meshData, z_knoppel_geomcentral, 4, 0.015);
    colormap(cbrewer('RdBu', 5));
    exportgraphics(gca, knoppel_geomcentral_image_name, Resolution=300);
end

end