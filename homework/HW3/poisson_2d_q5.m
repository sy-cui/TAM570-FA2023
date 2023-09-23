%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Poisson 2D Q5
rhs_func = @(x, y) 0*x + 1;
bc = 'dddd';
ns = 2.^(1:9);
ns = 2:2:64;
n2s = 2.*ns;
max_vals = ns;
for i = 1:length(ns);
    N = ns(i);
    N2 = n2s(i);
    disp(sprintf("Solving 2D N = %d", N))
    [~, ~, ub] = poisson2d(N, N, bc, rhs_func);
    [~, ~, u2b] = poisson2d(N2, N2, bc, rhs_func);
    max_diffs(i) = abs(max(max(abs(ub))) - max(max(abs(u2b))));
end;

figure('Units', 'inches', 'Position', [0 0 4 4])
    semilogy(ns, max_diffs, '-k', lw, 2)
    xlabel('$N$', intp, ltx); 
    ylabel('$L_\infty$ error', intp, ltx);
    pos = get(gcf, 'Position');
    set(gca, fs, 12, lw, 1, 'fontname', 'serif')
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, '2d_05_contour.pdf', '-dpdf')

