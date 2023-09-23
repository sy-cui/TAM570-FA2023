%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

%% Poisson 2D testing

% Q1
rhs_func{1} = @(x, y) 0*x + 1;
soln_func{1} = @(x, y) -0.5 * (x.^2 - 1);
bc{1} = 'ddnn';

% Q2
rhs_func{2} = @(x, y) 0*x + 1;
soln_func{2} = @(x, y) -0.5 * (y.^2 - 1);
bc{2} = 'nndd';

% Q3
rhs_func{3} = @(x, y) 0.5*pi^2*cos(0.5*pi*x).*cos(0.5*pi*y);
soln_func{3} = @(x, y) cos(0.5*pi*x).*cos(0.5*pi*y);
bc{3} = 'dddd';

% Q4
rhs_func{4} = @(x, y) 1.25*pi^2*cos(pi*x).*cos(0.5*pi*y);
soln_func{4} = @(x, y) cos(pi*x).*cos(0.5*pi*y);
bc{4} = 'nndd';

ns = [2 4 8 16 32 48 64 84 96 112 128];
errors = ns;

Q = 1
for i = 1:length(ns);
    N = ns(i);
    disp(sprintf("Solving 2D N = %d", N))
    [x_grid, y_grid, ub] = poisson2d(N, N, bc{Q}, rhs_func{Q});
    soln = soln_func{Q}(x_grid, y_grid);
    
    hx = x_grid(2, 1) - x_grid(1, 1);
    hy = y_grid(1, 2) - y_grid(1, 1);

    % errors(i) = sqrt(sum(sum(abs(soln - ub))) * hx * hy); % l_1
    % errors(i) = sqrt(sum(sum((soln - ub).^2)) * hx * hy); % l_2
    errors(i) = max(max(abs(soln - ub))); % l_inf
end;

figure('Units', 'inches', 'Position', [0 0 8 4])
    subplot(1, 2, 2)
    semilogy(ns, errors, '-k', lw, 2);
    xlabel('$N$', intp, ltx); 
    ylabel('$L_\infty$ error', intp, ltx);

    subplot(1, 2, 1)
    n_plot = 32;
    [x_grid, y_grid, ub] = poisson2d(n_plot, n_plot, bc{Q}, rhs_func{Q});
    surf(x_grid, y_grid, ub)
    pbaspect([1 1 1])
    colorbar
    xlabel('$x$', intp, ltx); 
    ylabel('$y$', intp, ltx);
    zlabel('$u(x, y)$', intp, ltx);
    pos = get(gcf, 'Position');
    set(gca, fs, 12, lw, 1, 'fontname', 'serif')
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, sprintf('2d_0%d_contour.pdf', Q), '-dpdf')
