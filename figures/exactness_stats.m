function [W1, W2, scaled_f_dist, mean_fiber_dist_cdf, fiber_dist_levels] = exactness_stats(meshData, z, f, S, degree, gij_sqrt)

N = size(f, 2);
[k, ik, ik_pinv] = fourier_ops(N, false, false);
f_sawtooth = ((N - 1) / (2 * pi)) * (ik_pinv .* (1 - abs(k) ./ (N / 2))).';
f_sawtooth = f_sawtooth .* exp(-ik.' .* angle(z));
f_sawtooth = ifft(f_sawtooth, [], 2, 'symmetric');

% normalize per fiber
f = f - mean(f, 2);
f = f ./ max(abs(f), [], 2);

f_sawtooth = f_sawtooth - mean(f_sawtooth, 2);
f_sawtooth = f_sawtooth ./ max(abs(f_sawtooth), [], 2);

f_dist = vecnorm(f - f_sawtooth, 2, 2);

scaled_f_dist = f_dist ./ vecnorm(f_sawtooth, 2, 2);

%% Distribution of density relative to extracted section
V2F = covariant_ops(meshData, degree);
z_colloc = V2F * z;
density = reshape(vecnorm(gij_sqrt.*S, 2, 1), 3 * meshData.nf, N);

angles = 2 * pi * ((0:(N - 1)) / N);
z_angle = mod(angle(z_colloc), 2 * pi);
fiber_dist = pi - abs(abs(mod(angles - z_angle + 2 * pi, 2 * pi)) - pi);

W1 = sum(density .* fiber_dist, 2) ./ sum(density, 2);
W2 = sum(density .* fiber_dist.^2, 2) ./ sum(density, 2);

fiber_dist_levels = pi * (0:(N - 1)) / N;
level_sets = fiber_dist < reshape(fiber_dist_levels, 1, 1, N);
fiber_dist_cdf = squeeze(mean(reshape(sum(density .* level_sets ./ sum(density, 2), 2), 3, meshData.nf, N), 1));

mean_fiber_dist_cdf = sum(meshData.areas .* fiber_dist_cdf, 1) / sum(meshData.areas, 1);

end