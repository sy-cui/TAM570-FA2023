% 08/30/2023
clc; clear all; close all;
format compact;

save_figure = true;

%% 1, 2. Plotting 

% Phase space definition
n = [10 20 40 80];
funcs = {@(x) sign(x), @(x) x, @(x) exp(cos(x)), @(x) cos(8*x+1)};
func_names = {'sign($x$)', '$x$', '$e^{\cos(x)}$', '$\cos(8x+1)$'};

l_inf = zeros(length(n), length(funcs));

% Subplot properties
subplot_i = floor(sqrt(length(n)));
subplot_j = ceil(length(n) / subplot_i);

for i = 1:length(n);
    N = n(i);
    M = 10 * N  + 10;

    [F, S, x, y, k] = compute_projection_matrices(N, M);

    for j = 1:length(funcs);
        f = funcs{j};
        u = f(x);           % Sampled points 
        u_k = F * u';       % Fourier coefficients
        u_s = (S * u_k)';   % Synthesized points
        u_a = f(y);         % Analytical function

        l_inf(i, j) = max(abs(u_a - u_s));

        figure(j, 'Units', 'inches', 'Position', [0 0 8 6]); box on;
        subplot(subplot_i, subplot_j, i); box on;
        hold on
        scatter(x, u, '-ob', 'linewidth', 1, 'DisplayName', 'Sampled points');
        plot(y, f(y), '-r', 'linewidth', 2, 'DisplayName', 'Original function');
        plot(y, u_s, '-k', 'linewidth', 1, 'DisplayName', 'Fourier interpolation');
        hold off
        title(sprintf('$N = %d$', N), 'interpreter', 'latex')
        xlim([-pi, pi])
        xlabel('$x$', 'interpreter', 'latex')
        ylabel('$f(x)$', 'interpreter', 'latex')
        set(gca, 'fontname', 'serif', 'fontsize', 12, 'linewidth', 1)
        legend('location', 'northwest', 'fontname', 'serif', 'fontsize', 12)
    end;
end;

if save_figure;
    for i = 1:length(funcs);
        figure(i);
        pos = get(gcf, 'Position');
        set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3) pos(4)]);
        print(gcf, sprintf('fourier_interpolation_%d.pdf', i), '-dpdf')
    end;
end;

pause
close all;

%% 3. Maximum point-wise error
l_inf
figure(length(funcs) + 1, 'Units', 'inches', 'Position', [0 0 5 4]); box on;
hold on
for j = 1:length(funcs);
    loglog(n, l_inf(:, j), '-o', 'linewidth', 2, 'DisplayName', func_names{j})
end;
hold off
xlabel('$N$', 'interpreter', 'latex')
ylabel('$L_\infty$ error', 'interpreter', 'latex')
legend('interpreter', 'latex', 'location', 'east')
set(gca, 'fontname', 'serif', 'fontsize', 12, 'linewidth', 1)

if save_figure;
    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3) pos(4)]);
    print(gcf, 'max_pointwise_error.pdf', '-dpdf')
end;

pause
close all

%% 4. Modal coefficients
N = 80;
[F, S, x, y, k] = compute_projection_matrices(N, 10 * N + 10);
indices = [1:N/2+1];
k_half = abs(k(indices));

figure(length(funcs) + 2, 'Units', 'inches', 'Position', [0 0 5 4]); box on;
hold on
for i = 1:length(funcs);
    f = funcs{i};
    u = f(x);           % Sampled points 
    u_k = F * u';       % Fourier coefficients

    abs_u_k_half = abs(u_k(indices));
    loglog(k_half, abs_u_k_half, 'linewidth', 2, 'DisplayName', func_names{i});
end;
hold off
xlabel('Fourier modes $|k|$', 'interpreter', 'latex')
ylabel('Fourier coefficients $|\hat{u}_k|$', 'interpreter', 'latex')
legend('interpreter', 'latex', 'location', 'east')
set(gca, 'fontname', 'serif', 'fontsize', 12, 'linewidth', 1)

if save_figure;
    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3) pos(4)]);
    print(gcf, 'modal_coefficients.pdf', '-dpdf')
end;
pause
close all



