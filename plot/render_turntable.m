function render_turntable(obj, filename, fps, duration, light)

if nargin < 5
    light = false;
end

set(obj, 'color', 'white');

nFrames = duration * fps;

if light
    h = camlight('left');
    lighting gouraud;
end

for j = 1:nFrames
    camorbit(360/nFrames, 0);
    if light
        camlight(h, 'left');
    end
    fr = getframe(obj);
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