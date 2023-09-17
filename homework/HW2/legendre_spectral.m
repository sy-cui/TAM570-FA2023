premeable;
%% 1. Interpolation
funcs = {@(x) exp(x), @(x) sign(x), @(x) 1 ./ (1 + 25 .* x.^2)};
func_names = {"$e^x$", "$sign(x)$", "$\frac{1}{1 + 25x^2}$"};

ny = 128;
n = 5:5:60;

h_inv = 0.5 * (ny - 1);
y = linspace(-1, 1, ny);                            % Interpolation points
interp = zeros(length(funcs), 2, length(n), ny);    % Interpolant buffer
l2_error = zeros(length(funcs), 2, length(n));      % L2 error norm buffer
soln = zeros(length(funcs), ny);                    % Original function buffer

for i = 1:length(funcs);
    soln(i, :) = funcs{i}(y);   % Compute analytical solution
end;

for i = 1:length(n);
    N = n(i);
    [zu, wu] = zwuni(N);
    [zg, wg] = zwgll(N);
    Ju = interp_mat(y, zu);
    Jg = interp_mat(y, zg);

    for j = 1:length(funcs); 
        interp(j, 1, i, :) = Ju * funcs{j}(zu);
        interp(j, 2, i, :) = Jg * funcs{j}(zg);
        l2_error(j, 1, i) = sqrt(sum((soln(j, :)' - squeeze(interp(j, 1, i, :))).^2) * h_inv);
        l2_error(j, 2, i) = sqrt(sum((soln(j, :)' - squeeze(interp(j, 2, i, :))).^2) * h_inv);
    end;

end;

for i = 1:length(funcs);
    figure(i, 'Units', 'inches', 'Position', [0 0 12 4]);
    subplot(1, 3, 1); box on;
        hold on;
        plot(y, soln(i, :), '-k', lw, 2, 'DisplayName', 'Analytical');
        for j = 1:2:length(n);
            plot(y, interp(i, 1, j, :), lw, 1, 'DisplayName', sprintf('N = %d', n(j)));
            % [z, ~] = zwuni(n(j));
            % scatter(z, funcs{i}(z), 'HandleVisibility', 'off');
        end;
        xlabel('$x$', intp, ltx); ylabel('$f(x)$', intp, ltx); xlim([-1, 1]);
        legend(intp, ltx, 'location', 'northwest');
        set(gca, lw, 2, fs, 12);
        hold off;
    subplot(1, 3, 2); box on;
        hold on;
        plot(y, soln(i, :), '-k', lw, 2, 'DisplayName', 'Analytical');
        for j = 1:2:length(n);
            plot(y, interp(i, 2, j, :), lw, 1, 'DisplayName', sprintf('N = %d', n(j)));
            % [z, ~] = zwgll(n(j));
            % scatter(z, funcs{i}(z), 'HandleVisibility', 'off');
        end;
        xlabel('$x$', intp, ltx); ylabel('$f(x)$', intp, ltx); xlim([-1, 1]);
        legend(intp, ltx, 'location', 'northwest');
        set(gca, lw, 2, fs, 12);
        hold off;
    subplot(1, 3, 3); grid on; box on;
        hold on;
        semilogy(n, l2_error(i, 1, :), '-k', lw, 2, 'DisplayName', 'Uniform')
        semilogy(n, l2_error(i, 2, :), '-r', lw, 2, 'DisplayName', 'GLL')
        xlabel('$n$', intp, ltx); ylabel('L2 error', intp, ltx);
        legend(intp, ltx, 'location', 'northwest');
        set(gca, lw, 2, fs, 12);
        hold off;
    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, sprintf('01_0%d_interp.pdf', i), '-dpdf')
end;
pause

%% 2. Legendre weighted residual method
ny = 256;
y = linspace(-1, 1, ny);
n = 2:24;
errors = n * 0.0;

for i = 1:length(n);
    N = n(i);
    [Ah,Bh,Ch,Dh,z,w] =  semhat(N);
    J = interp_mat(y, z);

    soln = sin(pi*(z+1)/4);
    f = ((pi/4)^2 + 1) * sin(pi*(z+1)/4);
    R = speye(N + 1); R = R(2:end, :);
    rhs = R * Bh * f;
    lhs = R * (Ah + Bh) * R';
    u = R' * (lhs \ rhs);

    errors(i) = max(abs(J * (u - soln)));
end;

figure(4, 'Units', 'inches', 'Position', [0 0 6 6]); box on; grid on;
    semilogy(n, errors, '-k', lw, 2)
    xlabel('$N$', intp, ltx); 
    ylabel('$L_\infty$ error', intp, ltx); 
    set(gca, lw, 1, fs, 14);
    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, '02_01_linf_convergence.pdf', '-dpdf')
pause

%% 4. Advetion example
c = 1; L = 1; 
nus = [0.01 0.005 0.001];
n = 5:5:150;
ny = 256;
y = linspace(0, L, ny);
errors = zeros(length(nus), length(n));

h = zeros(1, length(n));

for i = 1:length(n);
    N = n(i);

    R = speye(N + 1); R = R(2:end-1, :);
    [z, w] = zwgll(N);
    z = 0.5 * (z + 1) * L;
    w = 0.5 * w * L;
    h(i) = min(z(2:end) - z(1:end-1));
    Dh = new_deriv_mat(z);
    Bh = diag(w);
    Ah = Dh' * Bh * Dh;
    Ch = Bh * Dh;
    rhs = R * w;
    J = interp_mat(y, z);

    % For analytical solution
    xi = y / L;

    for j = 1:length(nus);
        % Compute spectral solution
        nu = nus(j);
        lhs = R * (nu * Ah + c * Ch) * R';
        u = R' * (lhs \ rhs);
        uy = J * u;

        % Compute analytical solution
        Pe = c * L / nu;
        us = L / c * (xi - (exp(Pe*(xi - 1)) - exp(-Pe))/(1 - exp(-Pe)));
        errors(j, i) = max(abs(reshape(us, ny, 1) - reshape(uy, ny, 1)));
    end;
end;

% Grid Peclet number
Pe_g = c * h ./ nus';
colors = {'black', 'red', 'green'};

yline = @(yval, varargin) plot(xlim, [yval yval], varargin{:});

figure(5, 'Units', 'inches', 'Position', [0 0 6 4]); box on; grid on;
    hold on;
    for i = 1:length(nus);
        semilogy(n, errors(i, :), '-', 'Color', colors{i}, 
            'linewidth', 2, 'DisplayName', sprintf('Error, $\\nu = %d$', nus(i)));
    end;
    for i = 1:length(nus);
        semilogy(n, Pe_g(i, :), '--', 'Color', colors{i}, 
            'linewidth', 2, 'DisplayName', sprintf('$Pe_g, \\nu = %d$', nus(i)));
    end;
    yline(2, '-.k')
    hold off;
    xlabel('$N$', intp, ltx); 
    ylabel('$L_\infty$ error', intp, ltx); 
    xlim([5, 150])
    legend('location', 'southwest', intp, ltx)
    set(gca, lw, 1, fs, 14);
    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
    print(gcf, '04_01_error_peg.pdf', '-dpdf')
