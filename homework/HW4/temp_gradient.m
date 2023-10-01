%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

ny = 128;
nx = 4000;
xm = 100;
ux_func = @(y) 1 - y.^2;
alpha = 0.01;
method = 'bdf2';
[x, y, soln, ~, ~, ~] = graetz(nx, ny, xm, ux_func, alpha, method);

yb = [-1 0 1]';
J = interp_mat(yb, y);
Dh = new_deriv_mat(y);
Gy = J * Dh * soln';

figure(1, 'Units', 'inches', 'Position', [2 2 4 4]); box on;
    hold on
    plot(x, Gy(1, :), '-k', lw, 1.5);
    plot(x, Gy(2, :), '-r', lw, 1.5);
    plot(x, Gy(3, :), '-b', lw, 1.5);
    hold off
    xlabel('$x$', intp, ltx); ylabel('$\partial T / \partial y$', intp, ltx);
    legend('$y = -1$', '$y = 0$', '$y = 1$', intp, ltx);
    set(gca, fs, 12, lw, 1.5, fn, 'serif');
    savefig_pdf('temp_gradient');
