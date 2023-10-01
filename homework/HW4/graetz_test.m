%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Definitions
ny = 128;
nx = 1000;
xm = 100;
ux_func = @(y) 1 - y.^2;
alpha = 0.01;
method = 'bdf2';
 
% Solve
[x, y, soln, lm, C, w] = graetz(nx, ny, xm, ux_func, alpha, method);

% Sampling
[y_uni, ~] = zwuni(64);
x_samp_idx = 1:50:nx;
x_samp = x(x_samp_idx);
J = interp_mat(y_uni, y);
[xx_samp, yy_samp] = ndgrid(x_samp, y_uni);
soln_samp = soln(x_samp_idx, :) * J';
eig_mode = C * exp(-lm*xx_samp).*(J*w)';

% Plotting
figure(1, 'Units', 'inches', 'Position', [2 2 4 4]); box on;
    mesh(xx_samp, yy_samp, soln_samp, lw, 1);
    xlim([0 xm]); ylim([-1 1]); zlim([0 1.25]);
    xlabel('$x$', intp, ltx); 
    ylabel('$y$', intp, ltx); 
    zlabel('$T(x, y)$', intp, ltx);
    view(40, 15);
    set(gca, lw, 1, fs, 12, fn, 'serif')
    savefig_pdf('graetz_mesh')

figure(2, 'Units', 'inches', 'Position', [2 2 8 4]); box on;
    subplot(1, 2, 1); mesh(xx_samp, yy_samp, eig_mode, lw, 1);
    title('Slowest-decaying eigenmode');
    subplot(1, 2, 2); mesh(xx_samp, yy_samp, abs(soln_samp-eig_mode), lw, 1);
    title('Absolute difference')
    for i=1:2;
        subplot(1, 2, i);
        xlim([0 xm]); ylim([-1 1]); zlim([0 1.25]);
        xlabel('$x$', intp, ltx); ylabel('$y$', intp, ltx); zlabel('$T(x, y)$', intp, ltx);
        view(40, 15);
        set(gca, lw, 1, fs, 12, fn, 'serif')
    end;
    savefig_pdf('graetz_decay')
