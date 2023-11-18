addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

function [y] = legendre_poly(x, N);
    % Evaluate Legendre polynomial of degree N at points x
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    if N == 0; y = x*0 + 1; return; end;

    [z, ~] = zwgl(N); z = z';
    x = reshape(x, length(x), 1);
    % Nth Lagrangian polynomial on GL points
    y = prod((x - z) ./ (1 - z), 2);

end;

%% 1. Plot Legendre polynomials and the new basis
x = linspace(-1, 1, 128);
kmax = 4;
colors = {'k','r','b','g','m'};

figure(1, 'Units', 'inches', 'Position', [2 2 8 4]);
    subplot(1,2,1);  box on; hold on;
    xlim([-1, 1]); ylim([-1 kmax+1]);% axis equal
    for k=0:kmax;
        yline(k, '--', 'color', colors{k+1}, lw, 1.5, 'HandleVisibility','off');
        plot(x, legendre_poly(x,k)+k, '-', 'color', colors{k+1}, 
            lw, 1.5, 'DisplayName', sprintf('$k=%d$',k));
    end;
    hold off;
    xlabel('$r$', intp, ltx); ylabel('$L_k(r) + k$', intp, ltx);
    set(gca, fs, 12, fn, 'serif', lw, 1.5);

    subplot(1,2,2);  box on; hold on;
    N = 10;
    xlim([-1, 1]); ylim([-1 kmax+1]);% axis equal
    for k=0:kmax;
        if k == 0;
            phi = 0.5 * (1 - x);
        elseif k == N;
            phi = 0.5 * (1 + x);
        else;
            phi = legendre_poly(x,k+1) - legendre_poly(x,k-1) + k;
        end;
        yline(k, '--', 'color', colors{k+1}, lw, 1.5, 'HandleVisibility','off');
        plot(x, phi, '-', 'color', colors{k+1}, 
            lw, 1.5, 'DisplayName', sprintf('$k=%d$',k));
    end;
    hold off;
    xlabel('$r$', intp, ltx); ylabel('$\phi_k(r) + k$', intp, ltx);
    set(gca, fs, 12, fn, 'serif', lw, 1.5);

savefig_pdf('basis_poly')

%% 2. Band structure of mass matrix
function [b] = phi_inner_prod(N1, N2);
    % Need quadrature to integrate exactly poly degree N1+N2+2
    [z, w] = zwgl(max(N1, N2)+3);
    b = sum(w.*(
        legendre_poly(z, N1+1).*legendre_poly(z, N2+1)
        + legendre_poly(z, N1-1).*legendre_poly(z, N2+1)
        + legendre_poly(z, N1+1).*legendre_poly(z, N2-1)
        + legendre_poly(z, N1-1).*legendre_poly(z, N2-1)
    ));
end;
N = 9;
B = zeros(N);
for i = 1:N; for j=1:N; B(i,j)=phi_inner_prod(i,j); end; end;
mask = abs(B) > 1e-14;

figure(2, 'Units', 'inches', 'Position', [2 2 14 4]); 
    subplot(1, 3, 1)
    imagesc(mask);
    axis equal;
    axis off;

    subplot(1, 3, 2)
    imagesc(kron(mask, mask))
    axis equal;
    axis off;

    subplot(1, 3, 3)
    imagesc(kron(mask, kron(mask, mask)))
    axis equal;
    axis off;
