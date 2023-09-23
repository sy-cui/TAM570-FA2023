function [x_grid, y_grid, z_grid, ub, ts] = poisson3d(nx, ny, nz, bc, rhs_func);

[Ah_x, Bh_x, Ch_x, Dh_x, z_x, w_x] = semhat(nx);
[Ah_y, Bh_y, Ch_y, Dh_y, z_y, w_y] = semhat(ny);
[Ah_z, Bh_z, Ch_z, Dh_z, z_z, w_z] = semhat(nz);

[x_grid, y_grid, z_grid] = ndgrid(z_x, z_y, z_z);

% Restriction matrices
Rx = speye(nx + 1); Ry = speye(ny + 1); Rz = speye(nz + 1);

% Parse boundary conditions
if length(bc) ~= 6; error("Unrecognized boudnary conditions."); end;
if lower(bc(1)) == 'd'; Rx = Rx(2:end, :); end;
if lower(bc(2)) == 'd'; Rx = Rx(1:end-1, :); end;
if lower(bc(3)) == 'd'; Ry = Ry(2:end, :); end;
if lower(bc(4)) == 'd'; Ry = Ry(1:end-1, :); end;
if lower(bc(5)) == 'd'; Rz = Rz(2:end, :); end;
if lower(bc(6)) == 'd'; Rz = Rz(1:end-1, :); end;

Ax = Rx * Ah_x * Rx';
Ay = Ry * Ah_y * Ry';
Az = Rz * Ah_z * Rz';
Bx = Rx * Bh_x * Rx';
By = Ry * Bh_y * Ry';
Bz = Rz * Bh_z * Rz';

[Sx, Lamx] = gen_eig_decomp(Ax, Bx);
[Sy, Lamy] = gen_eig_decomp(Ay, By);
[Sz, Lamz] = gen_eig_decomp(Az, Bz);

lx = full(diag(Lamx)); ly = full(diag(Lamy)); lz = full(diag(Lamz));

eval_inv = 1 ./ (
    reshape(lx, [length(lx), 1, 1]) 
    + reshape(ly, [1, length(ly), 1]) 
    + reshape(lz, [1, 1, length(lz)]) 
);

rhs = rhs_func(x_grid, y_grid, z_grid); 
rhs = tensor3(Bh_z, Bh_y, Bh_x, rhs);
rhs = tensor3(Rz, Ry, Rx, rhs);

tic(); 
u = tensor3(Sz', Sy', Sx', rhs);
u = eval_inv .* u;
u = tensor3(Sz, Sy, Sx, u);
ts = toc();

% Prolongate u
ub = tensor3(Rz', Ry', Rx', u);