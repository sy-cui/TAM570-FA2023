%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

%% Poisson 2D testing
rhs_func = @(x, y) 2*pi*pi*sin(pi * x).*cos(pi * y);
soln_func = @(x, y) sin(pi * x).*cos(pi * y);
bc = 'ddnn';

ns = [3:2:64];
errors = ns;

for i = 1:length(ns);
    N = ns(i);
    [x_grid, y_grid, ub] = poisson2d(N, N, bc, rhs_func);

    soln = soln_func(x_grid, y_grid);
    
    hx = x_grid(2, 1) - x_grid(1, 1);
    hy = y_grid(1, 2) - y_grid(1, 1);

    % errors(i) = sqrt(sum(sum(abs(soln - ub))) * hx * hy); % l_1
    % errors(i) = sqrt(sum(sum((soln - ub).^2)) * hx * hy); % l_2
    errors(i) = max(max(abs(soln - ub))); % l_inf
    
end;
figure
loglog(ns, errors, '-k', lw, 2);
xlabel('$N$', intp, ltx); 
ylabel('$L_\infty$ error', intp, ltx);
set(gca, fs, 16, lw, 1)

% contourf(x_grid, y_grid, ub)
% pbaspect([1 1 1])
% colorbar