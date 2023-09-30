%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Definitions
ny = 128;
nx = 100;
xm = 100;
ux_func = @(y) 1 - y.^2;
alpha = 0.01;
method = 'bdf2';
 
% Solve
[x, y, soln, lm, C, w] = graetz(nx, ny, xm, ux_func, alpha, method);

% Sampling
[y_uni, ~] = zwuni(64);
x_samp_idx = 1:4:nx;
x_samp = x(x_samp_idx);
J = interp_mat(y_uni, y);
[xx_samp, yy_samp] = ndgrid(x_samp, y_uni);
soln_samp = soln(x_samp_idx, :) * J';
eig_mode = C * exp(-lm*xx_samp).*(J*w)';


% Plotting
figure(1, 'Units', 'inches', 'Position', [2 2 10 5]); box on;
    ax = zeros(1, 2);
    ax(1) = subplot(1, 2, 1); ax(2) = subplot(1, 2, 2);
    subplot(1, 2, 1); mesh(xx_samp, yy_samp, soln_samp, lw, 1);
    subplot(1, 2, 2); mesh(xx_samp, yy_samp, eig_mode, lw, 1);

    for i=1:2;
        subplot(1, 2, i);
        xlim([0 xm]); ylim([-1 1]); zlim([0 1.25]);
        xlabel('$x$', intp, ltx); ylabel('$y$', intp, ltx); zlabel('$T(x, y)$', intp, ltx);
        view(50, 20);
        set(gca, lw, 1, fs, 12)

    end;



