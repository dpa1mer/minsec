function [k, ik, ik_pinv, d, d2] = fourier_ops(N, discrete, use_conj_sym)

if nargin < 3
    use_conj_sym = false;
end

assert(round(N / 2) == N / 2);
k = [0:(N/2), (-N/2 + 1 : -1)].';
if use_conj_sym
    k = k(1:(N/2));
end

ik = 1i .* k;
ik_pinv = ([1; ik(2:end)]).^(-1);

if discrete
    d = (exp(2 * pi * ik / N) - 1) * (N / (2 * pi));
    d2 = 4 * sin(2 * pi * k / (2 * N)).^2 * (N / (2 * pi)).^2;
else
    d = ik;
    d2 = k.^2;
end

end