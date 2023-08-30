%% Finding roots of Bessel function of the first kind
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
plot(1:n, bessel_roots, 'o-k', 'linewidth', 1.5)
xlabel('Root number $k$', 'interpreter', 'latex')
ylabel('Number of iterations')
set(gca, 'fontsize', 16, 'fontname', 'TimesNewRoman')
set(gcf,
    'PaperUnits', 'inches',
    'PaperSize', [pos(3), pos(4)])

print(gcf, '02_01_n_Pe_comp.pdf', '-dpdf')

%% Quadrature rule
% n = 8;
% result = zeros(1, n);
% for k = 1:n
%     [z, w] = zwgll(2^k); 
%     z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
%     result(k) = sum(w .* besselj(0, bessel_roots(45) .* z).^2);
%     if k > 1
%         disp(abs(result(k) - result(k - 1)));
%     end
    
% end
% plot(2.^(1:n), result); pause

%% Plot result at z = 0
% ng = 100;   % Number of grid points spanning R
% N = 50;     % Maximum N
% r = linspace(0, 1, ng);
% [z, w] = zwgll(256); 
% z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
% T = zeros(N, ng);
% dT = zeros(N, ng);

% figure(2); hold on;
% figure(3); hold on;
% for k = 1:N;
%     j0 = besselj(0, bessel_roots(k) .*z);
%     sh_by_ch = tanh(bessel_roots(k) * 20);
%     beta_k = sum(w.*j0.*z) / sum(w.*j0.^2.*z);
%     if k > 1;
%         T(k, :) = T(k-1, :) + beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
%         dT(k, :) = T(k-1, :) - beta_k * besselj(0, bessel_roots(k)*r);
%     else;
%         T(k, :) = beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
%         dT(k, :) = -beta_k * besselj(0, bessel_roots(k)*r);
%     end;
%     % figure(2); plot(r, T(k, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))
%     % figure(3); plot(r, dT(k, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))    
% end;
% % figure(2); plot(r, T(end, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))
% % figure(3); plot(r, dT(end, :), 'linewidth', 1.5, 'DisplayName', strcat('N = ', num2str(k)))  
% % figure(2); legend();
% % figure(3); legend();
% save("temperature.mat", "T")

%% Convergence
% load("temperature.mat")
% M = max(T, [], 2)';
% err = abs(M(1:end - 1) - M(end));
% figure(4)
% hold on
% loglog(1:49, err)
% loglog(1:49, (1:49).^(-1.5) / 10)   % ~1.5 order convergence
% hold off
% pause

% L dependence
ng = 100;   % Number of grid points spanning R
N = 50;     % Maximum N
L = 20;     % Maximum L
ln = 100;
L = linspace(1, L, ln);
r = linspace(0, 1, ng);
[z, w] = zwgll(256); 
z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]
T = zeros(ln, ng);

for l = 1:ln;
    for k = 1:N;
        j0 = besselj(0, bessel_roots(k) .*z);
        sh_by_ch = tanh(bessel_roots(k) * L(l));
        beta_k = sum(w.*j0.*z) / sum(w.*j0.^2.*z);
        T(l, :) = T(l, :) + beta_k / bessel_roots(k) * besselj(0, bessel_roots(k)*r) * sh_by_ch;
    end;
end;
M50 = max(T, [], 2);
figure(5);
plot(L, M50)
pause
    