%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Definitions
ny = 128;
nx = 400;
xm = 200;
ux_func = @(y) 1 - y.^2;
alpha = 0.01;
method = 'bdf2';
 
% Solve
[x, y, soln, lm, C] = graetz(nx, ny, xm, ux_func, alpha, method);

% Sampling
[y_uni, ~] = zwuni(64);
x_samp_idx = 1:4:nx;
x_samp = x(x_samp_idx);
J = interp_mat(y_uni, y);
[xx_samp, yy_samp] = ndgrid(x_samp, y_uni);
soln_samp = soln(x_samp_idx, :) * J';


% Plotting
figure(1, 'Units', 'inches', 'Position', [2 2 5 5]); box on;
    mesh(xx_samp, yy_samp, soln_samp, lw, 1);
    xlim([0 xm]); ylim([-1 1]); zlim([0 1.1]);
    hold on
    mesh(xx_samp, yy_samp, C*exp(-lm * xx_samp))
    hold off




