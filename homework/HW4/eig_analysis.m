%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Eigen analysis of Graetz problem via Legendre approximation
ny = 256;
alpha = 0.01;

[Ah, Bh, ~, ~, y, ~] = semhat(ny);
R = speye(ny + 1)(2:end-1, :);
ux = 1 - y.^2;

A = alpha*R*Ah*R';
B = R*(Bh.*ux)*R'; 
[S, Lam] = gen_eig_decomp(A, B);

% Eigenfunction comparison
k = 1:4;
S_diff = sin(0.5*pi*(y+1)*k);

figure(1, 'Units', 'inches', 'Position', [2 2 10 5])
    subplot(1, 2, 1); hold on; box on;
    subplot(1, 2, 2); hold on; box on;
    for i=k;
        subplot(1, 2, 1);
        plot(y, S_diff(:,i), lw, 1)
        xlim([-1, 1]); ylim([-1.5, 1.5])
        xlabel('$y$', intp, ltx); 
        ylabel('$w(y)$', intp, ltx); 
        title('$u(y) = 1$', intp, ltx);
        subplot(1, 2, 2);
        plot(y, [0; -S(:, i); 0], lw, 1)
        xlim([-1, 1]); ylim([-1.5, 1.5])
        xlabel('$y$', intp, ltx); 
        ylabel('$w(y)$', intp, ltx); 
        title('$u(y) = 1 - y^2$', intp, ltx);
    end;
    subplot(1, 2, 1); hold off; set(gca, fs, 12, fn, 'serif', lw, 1.5);
    subplot(1, 2, 2); hold off; set(gca, fs, 12, fn, 'serif', lw, 1.5);
    savefig_pdf('eigenfunction_comp')

% Eigenvalue comparison
k = 1:8;
lam = full(diag(Lam))'(k);
lam_diff = alpha * 0.25 * pi^2 * k.^2;
figure(2, 'Units', 'inches', 'Position', [2 2 5 5])
    box on; hold on;
    plot(k, lam_diff, 'ok', lw, 1.5)
    plot(k, lam, 'or', lw, 1.5)
    plot(k, abs(lam - lam_diff), 'ob', lw, 1.5)
    hold off;
    legend('$u(y) = 1$', '$u(y) = 1-y^2$', 'Difference', intp, ltx, 'location', 'northwest')
    xlabel('Eigenmode $k$', intp, ltx);
    ylabel('Eigenvalue $\lambda_k$', intp, ltx)
    set(gca, fs, 12, fn, 'serif', lw, 1.5);
    savefig_pdf('eigenvalue_comp')
    
