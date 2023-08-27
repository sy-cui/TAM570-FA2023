%% Problem 1
figure(1)
hold on

h = 2.^(-(0:52));

loglog(h, h.^2, '-k', 'linewidth', 2, 'DisplayName', 'Quadratic truncation error');
legend()
xlabel("$h$", 'interpreter', 'latex')
ylabel("Error")
xticks(10.^(-(0:2:16)))
yticks(10.^(-(0:2:16)))
ylim([1e-16 1e16])


% Problem 2
x = linspace(0, pi / 2, 50)';
e_d1u = max(abs((sin(x + h) - sin(x - h)) ./ (2 * h) - cos(x)), [], 1);
e_d2u = max(abs((sin(x + h) + sin(x - h) - 2 * sin(x)) ./ h.^2 + sin(x)), [], 1);

figure(1, 'Units', 'inches')
loglog(h, e_d1u, '-r', 'linewidth', 2, 'DisplayName', 'First-derivative observed error')
loglog(h, e_d2u, '-b', 'linewidth', 2, 'DisplayName', 'Second-derivative observed error')

loglog(h, 1.1e-16 ./ h, '--k', 'linewidth', 2, 'HandleVisibility','off')
loglog(h, 4.4e-16 ./ h.^2, '--k', 'linewidth', 2, 'HandleVisibility','off')
yticks(10.^((-16:2:16)))

pos = get(gcf, 'Position');
set(gca, 
    'fontsize', 16, 
    'fontname', 'TimesNewRoman')
set(gcf,
    'PaperUnits', 'inches',
    'PaperSize', [pos(3), pos(4)])

% print(gcf, '01_01_round_off.pdf', '-dpdf')

pause
