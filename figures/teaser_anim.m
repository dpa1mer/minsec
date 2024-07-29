function teaser_anim(datadir, datafile, quadfile)

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

%% Plot everything and set alpha values

% Plot section volume
fig_vol = figure(Color='w', Position = [0 0 1500 1500]); ax_vol = gca;
patch(Faces=meshData.faces, Vertices=meshData.verts, FaceColor=[0.8 0.8 0.8], EdgeColor='none');
plot_integral_curves(meshData, z, 4);
view(3); axis image vis3d off;
hold on;
h_vol = plot_section_volume(meshData, [1; 1; r^(-1)] .* S, Scale=.2 / (2 * pi), Repeat=4, Shift=0);
% clim auto;
alim manual;
% colormap parula;
% freezeColors;
alpha_vol = h_vol.FaceVertexAlphaData;
h_vol.FaceVertexAlphaData = 0.0 * alpha_vol;
movegui(fig_vol, "center");


% Plot section surface
fig_surf = figure(Color='w', Position = [0 0 1500 1500]); ax_surf = gca;
h_bg = patch(Faces=meshData.faces, Vertices=meshData.verts, FaceColor=[0.8 0.8 0.8], EdgeColor='none');
view(3); axis image vis3d off;
hold on;

h_surf = plot_section_surface(meshData, z, 0, .2, 4);
% colormap hsv;
% freezeColors;

% Plot field values
alpha_surf = h_surf.FaceVertexAlphaData;
alim manual
% h_surf.FaceVertexAlphaData = 0 * alpha_surf;
[h_lift, h_quiv] = plot_section_annotations(meshData, z, [4000; 8000], 0, .2, 4);
% alpha_marker = h_lift.MarkerEdgeAlpha;
% h_lift.MarkerEdgeAlpha = 0;

% Plot integral curves
h_int = plot_integral_curves(meshData, z, 4);
h_int.EdgeAlpha = 0;

% Plot quad mesh
% fig_quad = figure(Color='w'); fig_quad.Position(3:4) = [1024 1024]; ax_quad = gca;
[V, F] = readOBJ(fullfile(datadir, quadfile), Quads=true);
h_quad = patch('Faces', F, 'Vertices', V - [0 0 0.005], 'FaceColor', [0.9 0.9 0.9]);
view(3); axis image vis3d off;
movegui(fig_surf, "center");


h_link = linkprop([ax_vol, ax_surf], {'CameraPosition', 'CameraUpVector'});

%% Do the animation
az_init = 0; el_init = 90;
view(ax_vol, az_init, el_init);

fps = 30;

vw_vol = VideoWriter(fullfile(datadir, "teaser_vol.mj2"), "Archival");
vw_vol.open();

vw_surf = VideoWriter(fullfile(datadir, "teaser_surf.mj2"), "Archival");
vw_surf.open();

% vw_quad = VideoWriter(fullfile(datadir, "teaser_quad.avi"), "Uncompressed AVI");
% vw_quad.open();

duration_in = 3.5;
duration_mid = 8;
duration_out = 3.5;
duration_total = 15;
n_frames_total = duration_total * fps;
d_az = 360 / n_frames_total;
el = [ease_in_out(90, 30, linspace(0, 1, round(duration_in * fps)).');
      30 * ones(round(duration_mid * fps), 1);
      ease_in_out(30, 90, linspace(0, 1, round(duration_out * fps)).')];
d_el = diff([90; el]);

vol_fade = [ease_in_out(0, 1, linspace(0, 1, round(duration_in * fps)).');
            ones(round(duration_mid * fps) + round(duration_out * fps), 1)];
surf_fade = [ones(round((duration_in + duration_mid/4) * fps), 1);
             ease_in_out(1, 0.1, linspace(0, 1, duration_mid / 4 * fps).');
             0.1 * ones(duration_mid / 2 * fps, 1);
             ease_in_out(0.1, 0, linspace(0, 1, duration_out * fps).')];
int_fade = [zeros((duration_in + duration_mid/2) * fps, 1);
            ease_in_out(0, 1, linspace(0, 1, duration_mid/4 * fps).');
            ones(duration_mid/4 * fps, 1);
            ease_in_out(1, 0, linspace(0, 1, duration_out * fps).')];
bg_fade = [ones((duration_in + duration_mid) * fps, 1);
           ease_in_out(1, 0, linspace(0, 1, duration_out * fps).')];
marker_fade = [ones((duration_in + duration_mid/2) * fps, 1);
               ease_in_out(1, 0, linspace(0, 1, duration_mid/4 * fps).');
               zeros((duration_mid/4 + duration_out) * fps, 1)];

for j = 1:n_frames_total
    if j > 1 && vol_fade(j) ~= vol_fade(j - 1)
        h_vol.FaceVertexAlphaData = vol_fade(j) * alpha_vol;
    end
    if j > 1 && surf_fade(j) ~= surf_fade(j - 1)
        h_surf.FaceVertexAlphaData = surf_fade(j) * alpha_surf;
    end
    if j > 1 && int_fade(j) ~= int_fade(j - 1)
        h_int.EdgeAlpha = int_fade(j);
    end
    if j > 1 && bg_fade(j) ~= bg_fade(j - 1)
        h_bg.FaceAlpha = bg_fade(j);
    end
    if j > 1 && marker_fade(j) ~= marker_fade(j - 1)
        h_lift.MarkerEdgeAlpha = marker_fade(j);
    end
    if j == ((duration_in + duration_mid * 3/4) * fps)
        h_quiv.Visible = "off";
    end
    camorbit(ax_vol, d_az, d_el(j));

    fr = getframe(fig_vol);
    vw_vol.writeVideo(fr);
    fr = getframe(fig_surf);
    vw_surf.writeVideo(fr);
    % fr = getframe(fig_quad);
    % vw_quad.writeVideo(fr);
end

vw_vol.close();
vw_surf.close();
% vw_quad.close();


function f = ease_in_out(f0, f1, t, t_max)
    if nargin < 4
        t_max = 1;
    end
    t = t / t_max;
    s = (1 - cos(pi * t)) / 2;
    f = (1 - s) * f0 + s * f1;
end

end