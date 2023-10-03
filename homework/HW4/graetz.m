function [x, y, soln, lam_min, coeff, evec] = graetz(nx, ny, xm, ux_func, alpha, method);

if ~strcmp(lower(method), 'bdf2') && ~strcmp(lower(method), 'cn');
    error('Unrecognized method.')
end;

x = linspace(0, xm, nx);
dx = xm / (nx - 1);

[Ah, Bh, ~, ~, y, ~] = semhat(ny);
R = speye(ny + 1)(2:end-1, :);      % Dirichlet on both sizes

[xx, yy] = ndgrid(x, y);
soln = 0*xx;                        % Solution buffer
soln(1, :) = 1;                     % Initial condition
ux = ux_func(y);                    % Convecting velocity

A = alpha*R*Ah*R';                  % L = -alpha d^2/dy^2
B = R*(Bh.*ux)*R';                  % x-vel integrated into mass matrix

[S, Lam] = gen_eig_decomp(A, B);
lam_min = full(Lam(1, 1));
coeff = ones(1, ny-1) * B * S(:, 1);
evec = [0; S(:, 1); 0];

if strcmp(lower(method), 'bdf2');   % BDF2
    A_bdf1 = dx*A + full(B);        % Neumann op for first step (BDF1)
    A_bdf2 = 2*dx*A + 3*full(B);    % Neumann op for BDF2

    for i = 2:nx;
        if i == 2;
            soln(i, 2:end-1) = A_bdf1 \ (B*soln(i-1, 2:end-1)');
        else;
            soln(i, 2:end-1) = A_bdf2 \ (B*(4*soln(i-1, 2:end-1)-soln(i-2, 2:end-1))');
        end;
    end;

else;                               % Crank-Nicolson
    Al_cn = 2*full(B) + dx*A;
    Ar_cn = 2*full(B) - dx*A; 

    for i = 2:nx;
        soln(i, 2:end-1) = Al_cn \ (Ar_cn*soln(i-1, 2:end-1)');
    end;
end;