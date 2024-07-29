function ts = render_surface(meshData, az, el, az_light, el_light, color_field, cmap)

if nargin < 7
    cmap = cbrewer('RdYlBu', 1024);
end

bg_color = [1 1 1];
fg_color = [.7 .75 .8];

fig = figure(Color=bg_color); fig.Position = [0 0 1024 1024];
if nargin > 5
    ts = tsurf(meshData.faces, meshData.verts, 'FaceColor', 'interp', 'FaceVertexCData', color_field, 'EdgeColor', 'none');
    colormap(cmap);
else
    ts = tsurf(meshData.faces, meshData.verts, 'FaceColor', fg_color, 'EdgeColor', 'none');
end
view(az, el); axis image vis3d off;

l = camlight(az_light, el_light);
set(fig, Color=bg_color)
set(ts, fsoft);
add_shadow(ts, l, Color=bg_color*0.8, BackgroundColor=bg_color, Fade='infinite');
apply_ambient_occlusion(ts,'AddLights',true,'SoftLighting',true);

end