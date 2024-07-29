function render_slices(filename, meshData, f, fps)

nv = meshData.nv;
nf = meshData.nf;

fig = figure;
if size(f, 1) == nv
    ptc = trisurf(meshData.tri, 'FaceVertexCData', f(:, 1), 'FaceColor', 'interp', 'EdgeColor', 'none');
elseif size(f, 1) == 3 * nf
    pltVerts = meshData.verts(meshData.faces.', :);
    pltFaces = reshape(1:3*nf, 3, nf).';
    ptc = patch('Faces', pltFaces, 'Vertices', pltVerts, ...
                'FaceVertexCData', f(:, 1), 'FaceColor', 'interp', 'EdgeColor', 'none');
end
view(2); axis image off; colormap viridis;

set(fig, 'color', 'white');

nFrames = size(f, 2);

for j = 1:nFrames
    ptc.FaceVertexCData = f(:, j);
    fr = getframe(fig);
    im = frame2im(fr);
    [imind,cm] = rgb2ind(im, 256);
    if j == 1
        imwrite(imind, cm, filename, 'LoopCount', Inf, 'DelayTime', 1/fps);
    else
        imwrite(imind, cm, filename, 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
%     pause(1/60);
end

end