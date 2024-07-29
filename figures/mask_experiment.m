function mask_experiment(outdir)

assert(isfolder(outdir));

[verts, faces] = load_mesh('head0-subdiv.off');
head = process_surface(verts, faces, Flat=false, NormalizeArea=1);

[S, C, f, z, admmData] = minsec(head, 64, .1, 2, .1);
render_surface(head, 0, 10, 10, 45); hold on;
plot_integral_curves(head, z, 2, 2000);
plot_singularities(head, z, 2, 0.02);
colormap(cbrewer('RdBu', 5));
save(fullfile(outdir, "head_unmasked.mat"), "head", "S", "C", "f", "z", "admmData");
exportgraphics(gca, fullfile(outdir, "head_unmasked.png"), Resolution=300);

%% Compute distances (heat method)
dist = dist_from_bdry(head);
dist_edges = 0.5 * (dist(head.edges(:, 1)) + dist(head.edges(:, 2)));

%% Hard mask
mask = find(dist_edges > .5);
render_surface(head, 0, 10, 10, 45, double(dist > .5));
exportgraphics(gca, fullfile(outdir, "head_mask.png"), Resolution=300);

[S, C, f, z, admmData] = minsec(head, 64, .1, 2, .1, ExcludedSingularities=mask);
render_surface(head, 0, 10, 10, 45); hold on;
plot_integral_curves(head, z, 2, 2000);
plot_singularities(head, z, 2, 0.02);
colormap(cbrewer('RdBu', 5));
save(fullfile(outdir, "head_masked.mat"), "head", "S", "C", "f", "z", "admmData", "dist", "dist_edges", "mask");
exportgraphics(gca, fullfile(outdir, "head_masked.png"), Resolution=300);

%% Soft mask
soft_mask = 0.5 * (tanh(30 * (dist - 0.5)) + 1);
render_surface(head, 0, 10, 10, 45, soft_mask);
exportgraphics(gca, fullfile(outdir, "head_softmask.png"), Resolution=300);
soft_mask = 0.5 * (soft_mask(head.edges(:, 1)) + soft_mask(head.edges(:, 2)));
lambda = 10 * soft_mask + 0.1;

[S, C, f, z, admmData] = minsec(head, 64, lambda, 2, .1);
render_surface(head, 0, 10, 10, 45); hold on;
plot_integral_curves(head, z, 2, 2000);
plot_singularities(head, z, 2, 0.02);
colormap(cbrewer('RdBu', 5));
save(fullfile(outdir, "head_softmasked.mat"), "head", "S", "C", "f", "z", "admmData", "dist", "soft_mask", "lambda");
exportgraphics(gca, fullfile(outdir, "head_softmasked.png"), Resolution=300);

end