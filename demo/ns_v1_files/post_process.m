addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

figure(1, 'Units', 'inches', 'Position', [2 2 6 4]);
    box on; hold on;
    load("with_curlcurl.mat");
    semilogy(ti, ee, '-k', lw, 1.5);
    load("without_curlcurl.mat");
    semilogy(ti, ee, '-r', lw, 1.5);
    hold off;
    legend("With curlcurl", "Without curlcurl", fn, 'serif');
    xlabel("Time"); ylabel("$||u - u_e\||_2$", intp, ltx)
    set(gca, fs, 16, fn, 'serif', lw, 1.5);
    savefig_pdf("kovasznay_error")

