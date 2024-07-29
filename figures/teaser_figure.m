function teaser_figure(datadir, datafile, quadfile)

[verts, faces] = readOBJ("models/clover.obj");
meshData = process_surface(verts, faces, Flat=true, NormalizeArea=1);

assert(isfolder(datadir));

datafile = fullfile(datadir, datafile);
if isfile(datafile)
    load(datafile, "meshData", "S", "C", "f", "z", "admmData", "r", "lambda");
else
    r = 0.1;
    lambda = 0.06;
    [S, C, f, z, admmData] = minsec(meshData, 64, lambda, 4, r, Visualize=true, Repeat=4, Scale=.02, Shift=0.5, Frequency=10, SubsampleRate=0.4);
    save(datafile, "meshData", "S", "C", "f", "z", "admmData", "r", "lambda");
end

%% Plot section volume
fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_vol = plot_section_volume(meshData, [1; 1; r^(-1)] .* S, Scale=.2 / (2 * pi), Repeat=4, Shift=0);
view(3); axis image vis3d off;
hold on; patch(Faces=meshData.faces, Vertices=meshData.verts, FaceColor=[0.8 0.8 0.8], EdgeColor='none');
exportgraphics(gca, fullfile(datadir, "teaser_volume.png"), Resolution=300);

%% Plot section surface
fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_surf = plot_section_surface(meshData, z, 0, .2, 4);
view(3); axis image vis3d off;
hold on; patch(Faces=meshData.faces, Vertices=meshData.verts, FaceColor=[0.8 0.8 0.8], EdgeColor='none');
exportgraphics(gca, fullfile(datadir, "teaser_surface.png"), Resolution=300);

%% Plot field values
fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_surf = plot_section_surface(meshData, z, 0, .2, 4);
view(3); axis image vis3d off;
alpha = h_surf.FaceVertexAlphaData;
alim manual
h_surf.FaceVertexAlphaData = 0.1*alpha;
hold on; patch(Faces=meshData.faces, Vertices=meshData.verts, FaceColor=[0.8 0.8 0.8], EdgeColor='none');
plot_section_annotations(meshData, z, [4000; 8000], 0, .2, 4);
exportgraphics(gca, fullfile(datadir, "teaser_annot.png"), Resolution=300);

%% Plot integral curves
fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_surf = plot_section_surface(meshData, z, 0, .2, 4);
view(3); axis image vis3d off;
alpha = h_surf.FaceVertexAlphaData;
alim manual
h_surf.FaceVertexAlphaData = 0.1*alpha;
hold on; h_int = plot_integral_curves(meshData, z, 4);
exportgraphics(gca, fullfile(datadir, "teaser_int.png"), Resolution=300);

%% Plot quad mesh
fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_surf = plot_section_surface(meshData, z, 0, .2, 4);
view(3); axis image vis3d off;
alpha = h_surf.FaceVertexAlphaData;
alim manual
h_surf.FaceVertexAlphaData = 0.1*alpha;
hold on; 
[V, F] = readOBJ(fullfile(datadir, quadfile), Quads=true);
% fig = figure(Color='w'); fig.Position(3:4) = [1024 1024];
h_mesh = patch('Faces', F, 'Vertices', V, 'FaceColor', [0.9 0.9 0.9]); view(3); axis image vis3d off;
exportgraphics(gca, fullfile(datadir, "teaser_quad_mesh.png"), Resolution=300);


end