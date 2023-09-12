%% 1. Finding roots of Bessel function of the first kind
n = 50;
g1 = 2; g2 = 2 + 0.1;
bessel_roots = zeros(1, n);
iters = zeros(1, n);
for k = 1:n
    r1 = g1; r2 = g2; tmp = r2; iter = 0;
    while abs(besselj(0, r2)) >= 1e-14
        tmp = r2;
        r2 = r2 - besselj(0, r2) * (r2 - r1) / (besselj(0, r2) - besselj(0, r1));
        r1 = tmp;
        iter = iter + 1;
    end 
    bessel_roots(k) = r1;
    iters(k) = iter;
    g1 = g1 + pi;
    g2 = g1 + 0.1;
end

figure(1)
x = 0:0.01:200;
plot(1:n, iters, 'o-k', 'linewidth', 1.5)
xlabel('Root number $k$', 'interpreter', 'latex')
ylabel('Number of iterations')
set(gca, 'fontsize', 12, 'fontname', 'TimesNewRoman')
set(gcf,'PaperUnits', 'inches','PaperSize', [3, 4])
pause

%% 2. Quadrature rule
n = 8;
k_test = 50;
result = zeros(1, n);
q_diff = zeros(1, n - 1);
for k = 1:n
    [z, w] = zwgll(2^k); 
    z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
    result(k) = sum(w .* besselj(0, bessel_roots(k_test) .* z).^2 .* z);
    if k > 1
        q_diff(k - 1) = abs(result(k) - result(k - 1));
    end
end
figure(2, 'Units', 'inches', 'Position', [0 0 9 4])
subplot(1, 2, 1)
plot(2.^(1:n), result,'-k', 'linewidth', 2);
xlabel('Number of Quadrature Points $(n)$', 'interpreter', 'latex')
ylabel('Value of Integral $(\int_0^{1} J_{0}(\xi_k x)^2 \, x \, dx)$', 'interpreter', 'latex')
set(gca, 'fontsize', 12, 'fontname', 'serif')

subplot(1, 2, 2)
semilogy(2.^(2:n), q_diff, '--k', 'linewidth', 2)
xlabel('Number of Quadrature Points $(n)$', 'interpreter', 'latex')
ylabel('Absolute difference from $n - 1$', 'interpreter', 'latex')
set(gca, 'fontsize', 12, 'fontname', 'serif')

pos = get(gcf, 'Position');
set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3) pos(4)])

% print(gcf, '03_02_quad_convergence.pdf', '-dpdf')

%% 3. Plot result at z = 0
ng = 128;   % Number of grid points spanning R
N = 50;     % Maximum N
r = linspace(0, 1, ng);
[z, w] = zwgll(256); 
z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
T = zeros(N, ng);
dT = zeros(N, ng);

figure(3, 'Units', 'inches', 'Position', [0 0 8 6]); hold on; box on;
figure(4, 'Units', 'inches', 'Position', [0 0 8 6]); hold on; box on;
for k = 1:N
    % j0 = besselj(0, bessel_roots(k) .*z);
    % sh_by_ch = tanh(bessel_roots(k) * 20);
    % beta_k = sum(w.*j0.*z) / sum(w.*j0.^2.*z);
    % if k > 1
    %     T(k, :) = T(k-1, :) + beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
    %     dT(k, :) = dT(k-1, :) - beta_k * besselj(0, bessel_roots(k)*r);
    % else
    %     T(k, :) = beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
    %     dT(k, :) = -beta_k * besselj(0, bessel_roots(k)*r);
    % end
    if k > 1
        [temp_T, temp_dT] = het(r, 0, k, 256, 20, bessel_roots);
        T(k, :) = T(k - 1, :) + temp_T;
        dT(k, :) = dT(k - 1, :) + temp_dT;
    else
        [T(k, :), dT(k, :)] = het(r, 0, k, 256, 20, bessel_roots);
    end
    if k < 10
        figure(3); plot(r, T(k, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))
        figure(4); plot(r, dT(k, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))    
    end
end
figure(3); plot(r, T(end, :), '--k', 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(N)))
figure(4); plot(r, dT(end, :), '--k', 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(N)))  

figure(3); 
legend('location',"northeastoutside"); 
xlabel('Radial position ($r$)', 'Interpreter',"latex"); 
ylabel('Temperature ($T$)', 'Interpreter',"latex")
pos = get(gcf, 'Position');
set(gca, 'fontsize', 14, 'fontname', 'serif')
set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
hold off
% print(gcf, '03_03_temp_profile.pdf', '-dpdf')

figure(4); legend('location',"northeastoutside"); 
xlabel('Radial position ($r$)', 'Interpreter',"latex"); 
ylabel('Temperature gradient ($\partial_z T$)', 'Interpreter', "latex")
set(gca, 'fontsize', 14, 'fontname', 'serif')
set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)])
hold off
% print(gcf, '03_04_dT_profile.pdf', '-dpdf')
% save("temperature.mat", "T")

%% 4. Convergence of maximum temperature
N = 50;
L = 20;
r = linspace(0, 1, 128); z = linspace(0, L, 128)';
T = 0 * z * r;
MT = zeros(1, N);
for k = 1:N;
    [temp_T, ~] = het(r, z, k, 256, L, bessel_roots);
    T = T + temp_T;
    MT(k) = max(max(T));
end;

save("max_temp.mat", "MT")

% load("temperature.mat")
% M = max(T, [], 2)';
% err = abs(M(1:end - 1) - M(end));

load("max_temp.mat")
err = abs(MT(1:end - 1) - MT(end));

figure(5)
hold on
loglog(1:49, err,'linewidth', 1.5, 'DisplayName','Actual Error')
loglog(1:49, (1:49).^(-1.5) / 10,'linewidth', 1.5, 'DisplayName',"Order 1.5")   % ~1.5 order convergence
legend()
%xlabel('N'); ylabel('|M_N - M_{50}|')
xlabel('Number of modes (N)', 'Interpreter',"latex"); ylabel('$|M_N - M_{50}|$', 'Interpreter',"latex")
set(gca, 'fontsize', 12)
hold off

%% 5. L dependence
ng = 100;   % Number of grid points spanning R
N = 50;     % Maximum N
L = 20;     % Maximum L
ln = 100;
L = linspace(0.1, L, ln);
r = linspace(0, 1, ng);
[z, w] = zwgll(256); 
z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
T = zeros(ln, ng);

for l = 1:ln
    for k = 1:N
        j0 = besselj(0, bessel_roots(k) .*z);
        sh_by_ch = tanh(bessel_roots(k) * L(l));
        beta_k = sum(w.*j0.*z) / sum(w.*j0.^2.*z);
        T(l, :) = T(l, :) + beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
    end
end
M50 = max(T, [], 2);
figure(6);
plot(L, M50,'linewidth', 1.5)
xlabel('Length (L)',"Interpreter","latex")
ylabel('$M_{50}$',"Interpreter","latex")
set(gca, 'fontsize', 12)
