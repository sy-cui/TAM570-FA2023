function [x, y, soln, et] = adv_diff_2d(
    nv,                     % [Nx Ny]: number of grid points
    nu,                     % Diffusivity
    fcx,                    % cx(x, y): x-advection speed function
    fcy,                    % cy(x, y): y-advection speed function
    f0,                     % u^0(x, y): Initial condition evaluation function 
    fsrc,                   % f(t, x, y): Source term function 
    Tend,                   % Total march time
    CFL,                    % Courant-Friedrich-Lewy number
    ht=0,                   % Time step size. This is only used if advection speed is zero. 
    bct='dddd',             % 'lrbt': Boundary condtions. Defaults to all Dirichlet
    bcv=[0 0 0 0],          % [l r b t]: Boundary values. Defaults to all homogeneous
    bp=[-1 1 -1 1],         % [l r b t]: Boundary locations. Defaults to \hat{\Omega}
    fa={false, @(t,x,y) 0}, % Analytical solution.
    tol=1e-14               % General tolerance. Defaults to 1e-14
    );

%% ==================== Input argument check ====================
validateattributes(nv, {'numeric'}, {'vector', 'numel', 2, 'integer', '>', 2});
validateattributes(nu, {'numeric'}, {'scalar', 'nonnegative'});
validateattributes(Tend, {'numeric'}, {'scalar', 'nonnegative'});
validateattributes(CFL, {'numeric'}, {'scalar', 'nonnegative'});
validateattributes(ht, {'numeric'}, {'scalar', 'nonnegative'});
validateattributes(bcv, {'numeric'}, {'vector', 'numel', 4});
validateattributes(bp, {'numeric'}, {'vector', 'numel', 4});
validateattributes(tol, {'numeric'}, {'scalar', 'positive'});

if ~ischar(bct);
    error("Input argument 'bct' must be of string type"); 
elseif length(bct) ~= 4;
    error("Input argument 'bct' must be of length 4"); 
elseif strcmp(lower(bct), 'nnnn');
    error("All Neumann boundary condition is currently not supported");
end;

for i=1:4; if ~strcmp(lower(bct(i)), {'n','d'}); 
    error(['Unrecognized boundary type: ', bct]);
end; end;
%% ==================== Input argument check ends ====================

%% ==================== Setup ====================
nx = nv(1); ny = nv(2);
Jac_x = 0.5 * (bp(2) - bp(1));  % Jacobian in x
Jac_y = 0.5 * (bp(4) - bp(3));  % Jacobian in y
[Ah_x, Bh_x, Ch_x, Dh_x, z_x, w_x] = semhat(nx);
[Ah_y, Bh_y, Ch_y, Dh_y, z_y, w_y] = semhat(ny);

[xh, yh] = ndgrid(z_x, z_y);
x = Jac_x * (xh + 1) + bp(1); 
y = Jac_y * (yh + 1) + bp(3); 
cx = fcx(x, y); cy = fcy(x, y);

use_analytical_soln = fa{1};
fsoln = fa{2};

% Compute time step
hx_min = min(z_x(2) - z_x(1));
hy_min = min(z_y(2) - z_y(1));
cx_max = max(max(abs(cx)));
cy_max = max(max(abs(cy)));
c_max = max(max(sqrt(cx.^2+cy.^2)));
diffusion = (nu > tol);     % Diffusion flag
advection = (c_max > tol);  % Advection flag

if ~diffusion && ~advection; 
    error('No diffusion or advection present');
end;

if ht < tol && ~advection;
    error('Advection velocity and ht cannot both be zero');
elseif ht < tol; % Use CFL
    ht = CFL * min(hx_min / cx_max, hy_min / cy_max);
else; % Use min between provided ht and CFL 
    ht = min(CFL * min(hx_min / cx_max, hy_min / cy_max), ht); 
end;
n_steps = ceil(Tend / ht);
ht = Tend / n_steps;

% Restriction matrices
bc_d_mask = char(num2cell(bct)) == 'd';
Rx = speye(nx + 1); Ry = speye(ny + 1);
if bc_d_mask(1); Rx = Rx(2:end, :); end;
if bc_d_mask(2); Rx = Rx(1:end-1, :); end;
if bc_d_mask(3); Ry = Ry(2:end, :); end;
if bc_d_mask(4); Ry = Ry(1:end-1, :); end;

% Operation matrices
Ah_x = Ah_x / Jac_x;
Ah_y = Ah_y / Jac_y;
Bh_x = Bh_x * Jac_x;
Bh_y = Bh_y * Jac_y;
Ax = full(Rx*Ah_x*Rx'); 
Ay = full(Ry*Ah_y*Ry'); 
Bx = Rx*Bh_x*Rx';
By = Ry*Bh_y*Ry';

% BDFk/EXTk coefficients
bdf1 = [1 -1];
bdf2 = [3 -4 1] / 2;
bdf3 = [11 -18 9 -2] / 6;
ext1 = [1];
ext2 = [2 -1];
ext3 = [3 -3 1];

% Fast diagonalization eigen-decomposition
if diffusion;
    [Sx, Lx] = gen_eig_decomp(Ax, Bx);
    [Sy, Ly] = gen_eig_decomp(Ay, By);
    Lx = full(diag(Lx)); Ly = full(diag(Ly));
    Linv = 1 ./ (nu * ht * (Lx + Ly') + bdf1(1));
else;
    Linv = 1 ./ (bdf1(1)*full(diag(Bx)).*full(diag(By))');
end;

% Solution buffer
soln_buff = zeros(4, nx+1, ny+1); 
soln_buff(1,:,:) = f0(x, y); 

% Advection - Source buffer
ams_buff = zeros(4, nx+1, ny+1);
ams_buff(4,:,:) = -fsrc(0,x,y);
ams_buff(3,:,:) = ams_buff(4,:,:);
ams_buff(2,:,:) = ams_buff(4,:,:);
ams_buff(1,:,:) = (
    squeeze(ams_buff(4,:,:)) 
    + cx.*(Dh_x*squeeze(soln_buff(1,:,:))) 
    + cy.*(squeeze(soln_buff(1,:,:))*Dh_y')
);

% RHS buffer (Restricted)
rhs_buff = Linv * 0;

% Inhomogeneous boundary buffers
bc_d_val = reshape(bc_d_mask,4,1) .* reshape(bcv,4,1);
bc_n_val = reshape(1-bc_d_mask,4,1) .* reshape(bcv,4,1);

% Dirichlet
bc_d = sparse(nx+1,ny+1);
bc_d(1, :) = bc_d_val(1)*ones(1,ny+1);
bc_d(end, :) = bc_d_val(2)*ones(1,ny+1); 
bc_d(:, 1) = bc_d_val(3)*ones(nx+1,1); 
bc_d(:, end) = bc_d_val(4)*ones(nx+1,1); 
H_d = sparse(bdf1(1)*tensor2(Bh_y,Bh_x,bc_d)+nu*ht*(
    tensor2(Bh_y,Ah_x,bc_d)+tensor2(Ah_y,Bh_x,bc_d)
));

% Neumann
bc_n = sparse(nx+1,ny+1);
bc_n(1, :) = bc_n_val(1)*ones(1,ny+1)*Bh_y;
bc_n(end, :) = bc_n_val(2)*ones(1,ny+1)*Bh_y; 
bc_n(:, 1) = Bh_x*bc_n_val(3)*ones(nx+1,1); 
bc_n(:, end) = Bh_x*bc_n_val(4)*ones(nx+1,1); 
bc_n(:,:) = nu * ht * bc_n;

%% ==================== Setup ends ====================

% Time stepping loop
tic();
for k = 1:n_steps;
    % fprintf('Progress: %5.2f%% \r', k / n_steps * 100);
    kc = mod(k, 4) + 1;
    km1 = mod(k - 1, 4) + 1;
    km2 = mod(k - 2, 4) + 1;
    km3 = mod(k - 3, 4) + 1;

    if k > 2; % BDF3
        rhs_buff(:,:) = tensor2(Ry, Rx, tensor2(Bh_y, Bh_x, -squeeze(
            bdf3(2)*soln_buff(km1,:,:)
            + bdf3(3)*soln_buff(km2,:,:)
            + bdf3(4)*soln_buff(km3,:,:)
            + ht*ext3(1)*ams_buff(km1,:,:)
            + ht*ext3(2)*ams_buff(km2,:,:)
            + ht*ext3(3)*ams_buff(km3,:,:)
        )) + bc_n - H_d);

    elseif k == 2; % BDF2
        rhs_buff(:,:) = tensor2(Ry, Rx, tensor2(Bh_y, Bh_x, -squeeze(
            bdf2(2)*soln_buff(km1,:,:)
            + bdf2(3)*soln_buff(km2,:,:)
            + ht*ext2(1)*ams_buff(km1,:,:)
            + ht*ext2(2)*ams_buff(km2,:,:)
        )) + bc_n - H_d);
        
    else; % BDF1
        rhs_buff(:,:) = tensor2(Ry, Rx, tensor2(Bh_y, Bh_x, -squeeze(
            bdf1(2)*soln_buff(km1,:,:)
            + ht*ext1(1)*ams_buff(km1,:,:)
        )) + bc_n - H_d);

    end;

    if diffusion;
        % Fast Diagonalization Solve
        soln_buff(kc,:,:) = tensor2(Ry', Rx', 
            tensor2(Sy, Sx, Linv.*tensor2(Sy', Sx', rhs_buff))
        );
        % Update Linv
        if k == 1;
            Linv(:,:) = 1 ./ (nu * ht * (Lx + Ly') + bdf2(1));
        elseif k == 2;
            Linv(:,:) = 1 ./ (nu * ht * (Lx + Ly') + bdf3(1));
        end;
    else;
        % Direct inversion
        soln_buff(kc,:,:) = tensor2(Ry', Rx', Linv.*rhs_buff);
        % Update Linv
        if k == 1;
            Linv(:,:) = Linv * bdf1(1) / bdf2(1);
        elseif k == 2;
            Linv(:,:) = Linv * bdf2(1) / bdf3(1);
        end;
    end;
    soln_buff(kc,:,:) = squeeze(soln_buff(kc,:,:)) + bc_d;

    % Update advection - source term
    ams_buff(kc,:,:) = -fsrc(k*ht, x, y);
    if advection;
        ams_buff(kc,:,:) = (squeeze(ams_buff(kc,:,:)) 
            + cx.*(Dh_x*squeeze(soln_buff(kc,:,:))) 
            + cy.*(squeeze(soln_buff(kc,:,:))*Dh_y')
        );
    end;

    % Update H*u_b and force analytical input if applicable
    if k == 1;
        H_d = sparse(bdf2(1)*tensor2(Bh_y,Bh_x,bc_d)+nu*ht*(
            tensor2(Bh_y,Ah_x,bc_d)+tensor2(Ah_y,Bh_x,bc_d)
        ));
        if use_analytical_soln;
            soln_buff(kc,:,:) = fsoln(k*ht,x,y);
        end;
    elseif k == 2;
        H_d = sparse(bdf3(1)*tensor2(Bh_y,Bh_x,bc_d)+nu*ht*(
            tensor2(Bh_y,Ah_x,bc_d)+tensor2(Ah_y,Bh_x,bc_d)
        ));
        if use_analytical_soln;
            soln_buff(kc,:,:) = fsoln(k*ht,x,y);
        end;
    end;

    % if k==0 || k==n_steps|| mod(k, ceil(n_steps /20)) == 0;
    %     figure(1)
    %     mesh(x, y, squeeze(soln_buff(kc,:,:)))
    %     pbaspect([1 1 1])
    %     % zlim([0 1])
    %     pause(0.1)
    % end;
    
end; % Time loop
et = toc();
soln = squeeze(soln_buff(kc,:,:));
