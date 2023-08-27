clc; clear all; close all;
%% Comparison for several n and Pe
% Pes = [32, 64, 96];
% ns = [16, 32, 64];

% figure(1, 'Units', 'inches', 'Position', [0 0 10 10]);

% for m = 1:length(Pes);
%     for n = 1:length(ns);
%         subplot(length(Pes), length(ns), (m - 1) * length(ns) + n)
%         [x, u_sim, u_soln] = ad1d(ns(n), Pes(m), 1, 1);
%         plot(x, u_soln, '-k', 'linewidth', 2, 'DisplayName', 'Analytical solution');
%         hold on;
%         plot(x, u_sim, '-r', 'linewidth', 2, 'DisplayName', 'Centered difference');
%         hold off;
%         xlabel('$x$', 'interpreter', 'latex')
%         ylabel('$u(x)$', 'interpreter', 'latex')
%         title_string = sprintf('$n = %d,~Pe = %d,~Pe_g = %d$', ns(n), Pes(m), Pes(m) / ns(n));
%         title(title_string, 'interpreter', 'latex')
%         legend('location', 'northwest', 'fontsize', 14, 'fontname', 'TimesNewRoman')
%         set(gca, 'fontsize', 14, 'fontname', 'TimesNewRoman');
%     end;
% end;

% % pos = get(gcf, 'Position');
% % set(gcf,
% %     'PaperUnits', 'inches',
% %     'PaperSize', [pos(3), pos(4)])
% % print(gcf, '02_01_n_Pe_comp.pdf', '-dpdf')

%% Relative error vs n
Pes = [1 10 100 1000];
ns = 2.^(1:12);
errors = 0 * Pes' * ns;

figure(2, 'Units', 'inches', 'Position', [0 0 6 4]); hold on;
set(gca, 'fontsize', 14, 'fontname', 'Serif');
xlabel('$n$', 'interpreter', 'latex')
ylabel('$\| u - \tilde{u} \|_{\infty} / \| u \|_{\infty}$', 'interpreter', 'latex')
        
for m = 1:length(Pes);
    for n = 1:length(ns);
        [x, u_sim, u_soln] = ad1d(ns(n), Pes(m), 1, 1);
        errors(m, n) = max(abs(u_sim - u_soln)) / max(abs(u_soln));
    end;
    loglog(ns, errors(m, :), 'linewidth', 1.5, 'DisplayName', sprintf('$Pe = %d$', Pes(m)))
end;
legend('interpreter', 'latex')

loglog(ns, 1e-2 ./ ns.^2, '-.k', 'linewidth', 1.5, 'HandleVisibility', 'off')
xline = @(xval, varargin) plot([xval xval], ylim, varargin{:});
for Pe = Pes;
    if Pe ~= 1;
        xline(Pe / 2, '--k', 'HandleVisibility', 'off');
        t = text(
            Pe / 2 * 1.18, 1e3, sprintf('$n = %d$', Pe / 2), 
            'fontsize', 12, 
            'interpreter', 'latex',
            'rotation', 270);
    end;
end;

box on

pos = get(gcf, 'Position');
set(gcf,
    'PaperUnits', 'inches',
    'PaperSize', [pos(3), pos(4)])
print(gcf, '02_02_n_convergence.pdf', '-dpdf')
pause
