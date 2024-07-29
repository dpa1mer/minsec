function dataset_statistics(datadir)

data_files = dir(fullfile(datadir, "*.mat"));
n_data = length(data_files);
[~, base_name] = fileparts(string({data_files(:).name}.'));


global_stats = table(RowNames=base_name);
fiber_dist_cdf_tbl = table;

W2_cdf_tbl = table;
n_bins = 200;
W2_cdf_tbl.bin = linspace(0, 10, n_bins).';

for j = 1:n_data
    the_file = data_files(j);
    filename = fullfile(the_file.folder, the_file.name);

    load(filename, "meshData", "admmData", "S", "f", "z");
    % Calculate evaluation statistics
    [~, W2, ~, mean_fiber_dist_cdf, fiber_dist_levels] = exactness_stats(meshData, z, f, S, 4, [1;1;0.2^(-1)]);
    
    global_stats.nv(j) = meshData.nv;
    global_stats.nf(j) = meshData.nf;
    global_stats.time(j) = admmData(end).wallTime;
    global_stats.iters(j) = admmData(end).iter;

    fiber_dist_cdf_tbl.dist_level = fiber_dist_levels.';
    fiber_dist_cdf_tbl.(base_name(j)) = reshape(mean_fiber_dist_cdf, [], 1);

    W2_hist_weights = reshape(repmat(reshape(meshData.areas, 1, meshData.nf), 3, 1), 3 * meshData.nf, 1);
    W2_hist_weights = W2_hist_weights / (3 * sum(meshData.areas));
    W2_cdf_tbl.(base_name(j)) = cumsum(histwv(W2, W2_hist_weights, 0, 10, n_bins));

    pointwise_curvature = meshData.star0lump \ meshData.angleDefect;
    interior_curvature = pointwise_curvature(meshData.intIdx);
    M_int = meshData.star0lump(meshData.intIdx, meshData.intIdx);

    global_stats.mean_min_angle(j) = mean(min(meshData.faceVertAngles, [], 2), 1);
    global_stats.min_min_angle(j) = min(meshData.faceVertAngles, [], 'all');
    global_stats.max_angle_defect(j) = max(abs(meshData.angleDefect(meshData.intIdx)));
    global_stats.mean_abs_angle_defect(j) = mean(abs(meshData.angleDefect(meshData.intIdx)));
    global_stats.rms_angle_defect(j) = sqrt(mean(meshData.angleDefect(meshData.intIdx).^2));
    global_stats.max_abs_curvature(j) = max(abs(interior_curvature));
    global_stats.mean_abs_curvature(j) = sum(M_int * abs(interior_curvature) / sum(diag(M_int)));
    global_stats.rms_curvature(j) = sqrt(sum(M_int * (interior_curvature.^2)) / sum(diag(M_int)));
    global_stats.mean_tri_area(j) = mean(meshData.areas);
    global_stats.std_tri_area(j) = std(meshData.areas);
    global_stats.min_tri_area(j) = min(meshData.areas);
    global_stats.mean_W2(j) = sum(W2_hist_weights .* W2);
end

writetable(global_stats, fullfile(datadir, "global_stats.csv"), WriteRowNames=true);
writetable(fiber_dist_cdf_tbl, fullfile(datadir, "fiber_dist_cdf.csv"));
writetable(W2_cdf_tbl, fullfile(datadir, "W2_cdf.csv"));


end