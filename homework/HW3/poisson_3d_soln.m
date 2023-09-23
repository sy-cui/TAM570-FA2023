addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

rhs_func = @(x, y, z) (0.75*pi^2).*cos(0.5*pi* x).*cos(0.5*pi*y).*cos(0.5*pi*z);
soln_func = @(x, y, z) cos(0.5*pi*x).*cos(0.5*pi*y).*cos(0.5*pi*z);
bc = 'dddddd';

ns = [4 6 8 12 16 24 32 48 64 96 128 256 384 512 640];
errors = ns;
ts = ns;

for i = 1:length(ns);
    N = ns(i);
    disp(sprintf("N = %d", N))
    [x_grid, y_grid, z_grid, ub, t] = poisson3d(N,N,N,'dddddd',rhs_func);
    soln = soln_func(x_grid, y_grid, z_grid);

    ts(i) = t;
    
    hx = x_grid(2, 1, 1) - x_grid(1, 1, 1);
    hy = y_grid(1, 2, 1) - y_grid(1, 1, 1);
    hz = z_grid(1, 1, 2) - z_grid(1, 1, 1);

    % errors(i) = sqrt(sum(sum(sum(abs(soln - ub)))) * hx * hy * hz); % l_1
    % errors(i) = sqrt(sum(sum(sum((soln - ub).^2))) * hx * hy * hz); % l_2
    errors(i) = max(max(max(abs(soln - ub)))); % l_inf
    
end;

save("benchmark.mat", "ns", "ts", "errors")
    