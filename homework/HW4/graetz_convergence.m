%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

run_spatial_convergence = false;
run_temporal_convergence = false;

% Spatial convergence
if run_spatial_convergence;

    ny_arr = [4 8 16 32 64 96 128 192 256];
    nx = 4000;
    xm = 20;
    ux_func = @(y) 1 - y.^2;
    alpha = 0.01;
    method = 'bdf2';

    ny_samp = 128;
    solns_y = zeros(length(ny_arr), ny_samp+1);
    [y_uni, ~] = zwuni(ny_samp);

    for i = 1:length(ny_arr);
        ny = ny_arr(i);
        disp(['Solving ny = ', num2str(ny)]);
        [x,y,soln,~,~,~] = graetz(nx, ny, xm, ux_func, alpha, method);
        J = interp_mat(y_uni, y);
        solns_y(i, :) = soln(end, :)*J';
    end;

    save('convergence_y.mat', 'solns_y', 'ny_arr')

end;

% Temporal convergence
if run_temporal_convergence; 

    nx_arr = [1 5 10 20 30 50 80 100 200 500 1000 2000 4000];
    ny = 256;
    xm = 20;
    ux_func = @(y) 1 - y.^2;
    alpha = 0.01;
    method = 'bdf2';

    ny_samp = 128;
    solns_x = zeros(length(nx_arr), ny_samp+1);
    [y_uni, ~] = zwuni(ny_samp);

    for i = 1:length(nx_arr);
        nx = nx_arr(i);
        disp(['Solving nx = ', num2str(nx)]);
        [x,y,soln,~,~,~] = graetz(nx, ny, xm, ux_func, alpha, method);
        J = interp_mat(y_uni, y);
        solns_x(i, :) = soln(end, :)*J';
    end;

    save('convergence_x.mat', 'solns_x', 'nx_arr')

end;

% Plotting
load('convergence_y.mat')
load('convergence_x.mat')
figure(1, 'Units', 'inches', 'Position', [2 2 10 5]);
    subplot(1, 2, 1); box on;
    linf_error = max(abs(solns_y(1:end-1, :) - solns_y(end, :)), [], 2);
    semilogy(ny_arr(1:end-1), linf_error, '-ok', lw, 1.5);
    xlabel('$N_y$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    set(gca, lw, 1, fs, 12)

    subplot(1, 2, 2); box on;
    linf_error = max(abs(solns_x(1:end-1, :) - solns_x(end, :)), [], 2);
    hold on;
    loglog(nx_arr(1:end-1), linf_error, '-ok', lw, 1.5);
    loglog(nx_arr(1:end-1), nx_arr(1:end-1).^(-2), '--r', lw, 1.5)
    hold off;
    xlabel('$N_x$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    set(gca, lw, 1, fs, 12)

    savefig_pdf('graetz_convergece')