addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Post processing
load("benchmark.mat", "ns", "ts", "errors")

figure(1, 'Units', 'inches', 'Position', [0 0 4 4]); box on;
    semilogy(ns(1:end-3), errors(1:end-3), '-k', lw, 2);
    xlabel('$N$', intp, ltx); 
    ylabel('$L_\infty$ error', intp, ltx);
    pos = get(gcf, 'Position');
    set(gca, fs, 12, lw, 1, 'fontname', 'serif')
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, '3d_error.pdf', '-dpdf')

mflops = 12*ns.^4./ts*1e-9;
P1 = polyfit(log(ns), log(ts), 1)(1);
P2 = polyfit(log(ns), log(mflops), 1)(1);
figure(2, 'Units', 'inches', 'Position', [0 0 9 4]); box on;
    subplot(1, 2, 1); box on;
    hold on
    loglog(ns, ts, '-k', lw, 2);
    loglog(ns, ns.^P1 .*1e-6, '--r', lw, 2)
    hold off
    xlabel('$N$', intp, ltx); 
    ylabel('Wall time [s]', intp, ltx);
    set(gca, fs, 12, lw, 1, 'fontname', 'serif')

    subplot(1, 2, 2); box on;
    hold on
    loglog(ns, mflops, '-k', lw, 2);
    loglog(ns, 1e-2*ns.^P2, '--r', lw, 2);
    hold off
    xlabel('$N$', intp, ltx); 
    ylabel('MFlops', intp, ltx);
    set(gca, fs, 12, lw, 1, 'fontname', 'serif')

    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, '3d_timing.pdf', '-dpdf')