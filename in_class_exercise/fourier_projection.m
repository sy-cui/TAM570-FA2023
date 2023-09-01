% 08/30/2023

format compact; box on;

Ns = [10 20 40 80];
for j = 1:4;
subplot(2, 2, j)
N = Ns(j);
M = 10 * N + 10;
k = [-N/2:(N/2-1)]';

h = 2 * pi / N; h_f = 2 * pi / M;
x = -pi + h * [0:N-1]';
x_f = -pi + h_f * [0:M-1]';

i = sqrt(-1);
F = exp(-i * k * x') / N;

I = exp(i * x_f * k');

f = sign(x); f(1) = 0;
f = x;
% f = exp(cos(x));
% f = cos(x + 2);

interp = I * F * f;
hold on
plot(x, f, '-k', 'linewidth', 2, 'DisplayName', 'Original function');
plot(x_f, interp, '-r', 'linewidth', 2, 'DisplayName', 'Fourier interpolation');
hold off
title(sprintf('N = %d', N))
xlim([-pi, pi])
end;
pause

