function [x_grid, y_grid, ub] = poisson2d(nx, ny, bc, rhs_func);

[Ah_x, Bh_x, Ch_x, Dh_x, z_x, w_x] = semhat(nx);
[Ah_y, Bh_y, Ch_y, Dh_y, z_y, w_y] = semhat(ny);

[x_grid, y_grid] = ndgrid(z_x, z_y);

% Restriction matrices
Rx = speye(nx + 1); Ry = speye(ny + 1);

% Parse boundary conditions
if length(bc) ~= 4; error("Unrecognized boudnary conditions."); end;
if lower(bc(1)) == 'd'; Rx = Rx(2:end, :); end;     % Left boundary is Dirichlet 
if lower(bc(2)) == 'd'; Rx = Rx(1:end-1, :); end;   % Right boundary is Dirichlet
if lower(bc(3)) == 'd'; Ry = Ry(2:end, :); end;     % Bottom boundary is Dirichlet
if lower(bc(4)) == 'd'; Ry= Ry(1:end-1, :); end;    % Top boundary is Dirichlet

Ax = Rx * Ah_x * Rx';
Ay = Ry * Ah_y * Ry';
Bx = Rx * Bh_x * Rx';
By = Ry * Bh_y * Ry';

[Sx, Lamx] = gen_eig_decomp(Ax, Bx);
[Sy, Lamy] = gen_eig_decomp(Ay, By);

lx = full(diag(Lamx)); ly = full(diag(Lamy));

eval_inv = 1 ./ (
    reshape(lx, [length(lx), 1]) 
    + reshape(ly, [1, length(ly)]) 
);

rhs = rhs_func(x_grid, y_grid);
rhs = tensor2(Bh_y, Bh_x, rhs);
rhs = tensor2(Ry, Rx, rhs);

u = tensor2(Sy', Sx', rhs);
u = eval_inv .* u;
u = tensor2(Sy, Sx, u);

% Prolongate u
ub = tensor2(Ry', Rx', u);