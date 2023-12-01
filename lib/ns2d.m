function [x,y,soln,err,et] = ns2d(
    nv,                     % [Nx Ny]: number of grid points
    mv,                     % [Mx My]: number of over-integration points
    Re,                     % Reynolds number
    Pr,                     % Prantl number
    Tend,                   % Total march time
    CFL,                    % Courant-Friedrich-Lewy number
    Vscale,                 % Characteristic velocity scaled. Used to compute dt with CFL
    ht=0,                   % Override time step size
    bp=[-1 1 -1 1],         % Boundary points. Defaults to [-1 1 -1 1]
    % bc={                    % Boundary conditions. Defaults to all periodic
    %     'pppp', [0 0 0 0],  % u
    %     'pppp', [0 0 0 0],  % v
    %     'pppp', [0 0 0 0],  % p
    % },            
    f0={@(t,x,y) 0,         % Initial condition for U, V, and T
        @(t,x,y) 0,
        @(t,x,y) 0},
    fq={@(t,x,y) 0,         % Source terms for U, V, and T
        @(t,x,y) 0,
        @(t,x,y) 0},
    fe={false, @(t,x,y) 0,  % Analytical solution flag and function for U, V, and T 
        false, @(t,x,y) 0,
        false, @(t,x,y) 0},
    tol=1e-14               % General tolerance. Defaults to 1e-14
    );

%% ==================== Input argument check ====================
validateattributes(nv, {'numeric'}, {'vector', 'numel', 2, 'integer', '>', 1});
validateattributes(mv, {'numeric'}, {'vector', 'numel', 2, 'integer', '>', 1});
validateattributes(Re, {'numeric'}, {'scalar', 'positive'});
validateattributes(Pr, {'numeric'}, {'scalar', 'positive'});
validateattributes(Tend, {'numeric'}, {'scalar', 'positive'});
validateattributes(CFL, {'numeric'}, {'scalar', 'positive'});
validateattributes(Vscale, {'numeric'}, {'scalar', 'positive'});
validateattributes(ht, {'numeric'}, {'scalar', 'nonnegative'});
validateattributes(bp, {'numeric'}, {'vector', 'numel', 4});
validateattributes(tol, {'numeric'}, {'scalar', 'positive'});

% if ~ischar(bct);
%     error("Input argument 'bct' must be of string type"); 
% elseif length(bct) ~= 4;
%     error("Input argument 'bct' must be of length 4"); 
% elseif strcmp(lower(bct), 'nnnn');
%     error("All Neumann boundary condition is currently not supported");
% end;

% for i=1:4; if ~strcmp(lower(bct(i)), {'n','d'}); 
%     error(['Unrecognized boundary type: ', bct]);
% end; end;
%% ==================== Input argument check ends ====================

%% ==================== Setup ====================

% Interpret function inputs
u_ini = f0{1}; v_ini = f0{2}; t_ini = f0{3};
u_src = fq{1}; v_src = fq{2}; t_src = fq{3};
u_ext = fe{2}; v_ext = fe{4}; t_ext = fe{6};
exact_u = fe{1}; exact_v = fe{3}; exact_t = fe{5};

% Set up basic parameters
nx = nv(1); ny = nv(2);
Jac_x = 0.5 * (bp(2) - bp(1));  % Jacobian in x
Jac_y = 0.5 * (bp(4) - bp(3));  % Jacobian in y
[Ah_x, Bh_x, Ch_x, Dh_x, xi_x, rho_x] = semhat(nx);
[Ah_y, Bh_y, Ch_y, Dh_y, xi_y, rho_y] = semhat(ny);

xh = Jac_x * (xi_x + 1) + bp(1);
yh = Jac_y * (xi_y + 1) + bp(3);
[x, y] = ndgrid(xh, yh);

Pe = Re * Pr;   % Thermal Peclet number

% Set up dealiasing
mx = mv(1); my = mv(2);
[eta_x, w_x] = zwgl(mx);
[eta_y, w_y] = zwgl(my);
Jx = interp_mat(eta_x, xi_x);
Jy = interp_mat(eta_y, xi_y);

% For plotting
n_uni = 32;
[zu_x, ~] = zwuni(n_uni); [zu_y, ~] = zwuni(n_uni);
zu_x = Jac_x * (zu_x + 1) + bp(1);
zu_y = Jac_y * (zu_y + 1) + bp(3);
Ju_x = interp_mat(zu_x, xh); Ju_y = interp_mat(zu_y, yh); 
[xu, yu] = ndgrid(zu_x, zu_y);

% Compute time step
hx_min = min(xi_x(2)-xi_x(1)); hy_min = min(xi_y(2)-xi_y(1));
if ht < tol;
    ht = CFL * min(hx_min, hy_min) / Vscale;
end;
n_steps = ceil(Tend / ht);
ht = Tend / n_steps;

% BDFk/EXTk coefficients
bdf1 = [1 -1];
bdf2 = [3 -4 1] / 2;
bdf3 = [11 -18 9 -2] / 6;
ext1 = [1];
ext2 = [2 -1];
ext3 = [3 -3 1];

# Restriction / Periodic boundary matrices
Ru_x = speye(nx+1)(2:end-1,:); Ru_y = speye(ny+1)(2:end-1,:);
Rv_x = speye(nx+1)(2:end-1,:); Rv_y = speye(ny+1)(2:end-1,:);
Rp_x = speye(nx+1)(2:end-1,:); Rp_y = speye(ny+1);
Rt_x = speye(nx+1)(2:end-1,:); Rt_y = speye(ny+1)(2:end-1,:);
Ru_x = r_periodic(nx); Ru_y = r_periodic(ny);
Rv_x = r_periodic(nx); Rv_y = r_periodic(ny);
Rp_x = r_periodic(nx); Rp_y = r_periodic(ny);
Rt_x = r_periodic(nx); Rt_y = r_periodic(ny);

% Operation matrices
Ah_x = Ah_x / Jac_x; Ah_y = Ah_y / Jac_y;
Bh_x = Bh_x * Jac_x; Bh_y = Bh_y * Jac_y;
Dh_x = Dh_x / Jac_x; Dh_y = Dh_y / Jac_y;
Au_x = Ru_x*Ah_x*Ru_x'; Au_y = Ru_y*Ah_y*Ru_y';
Bu_x = Ru_x*Bh_x*Ru_x'; Bu_y = Ru_y*Bh_y*Ru_y';
Av_x = Rv_x*Ah_x*Rv_x'; Av_y = Rv_y*Ah_y*Rv_y';
Bv_x = Rv_x*Bh_x*Rv_x'; Bv_y = Rv_y*Bh_y*Rv_y';
Ap_x = Rp_x*Ah_x*Rp_x'; Ap_y = Rp_y*Ah_y*Rp_y';
Bp_x = Rp_x*Bh_x*Rp_x'; Bp_y = Rp_y*Bh_y*Rp_y';
At_x = Rt_x*Ah_x*Rt_x'; At_y = Rt_y*Ah_y*Rt_y';
Bt_x = Rt_x*Bh_x*Rt_x'; Bt_y = Rt_y*Bh_y*Rt_y';
Bh_xy = (Jac_y*Jac_x)*rho_x.*rho_y';  % F.*Bh_xy = tensor2(Bh_y, Bh_x, F)

JtBm_x=Jx'*diag(w_x)*Jac_x; JtBm_y=Jy'*diag(w_y)*Jac_y; 
JD_x = Jx*Dh_x; JD_y = Jy*Dh_y;

% Eigen decomposition
[Su_x,Lu_x]=gen_eig_decomp(Au_x,Bu_x); [Su_y,Lu_y]=gen_eig_decomp(Au_y,Bu_y);
[Sv_x,Lv_x]=gen_eig_decomp(Av_x,Bv_x); [Sv_y,Lv_y]=gen_eig_decomp(Av_y,Bv_y);
[Sp_x,Lp_x]=gen_eig_decomp(Ap_x,Bp_x); [Sp_y,Lp_y]=gen_eig_decomp(Ap_y,Bp_y);
[St_x,Lt_x]=gen_eig_decomp(At_x,Bt_x); [St_y,Lt_y]=gen_eig_decomp(At_y,Bt_y);
Su_x = Ru_x'*Su_x; Su_y = Ru_y'*Su_y;
Sv_x = Rv_x'*Sv_x; Sv_y = Rv_y'*Sv_y;
Sp_x = Rp_x'*Sp_x; Sp_y = Rp_y'*Sp_y;
St_x = Rt_x'*St_x; St_y = Rt_y'*St_y;
Lu_inv = 1.0 ./ (ht/Re*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf1(1));
Lv_inv = 1.0 ./ (ht/Re*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf1(1));
Lt_inv = 1.0 ./ (ht/Pe*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf1(1));
Lp_inv = 1.0 ./ (full(diag(Lp_x)) + full(diag(Lp_y))'); Lp_inv(1,1) = 0;

% Helper functions
function vorticity = compute_vorticity(velx, vely);
    vorticity = Dh_x*vely - velx*Dh_y';
end;

function [Curl_x, Curl_y] = curl(field_z);
    Curl_x = Bh_xy.*(field_z*Dh_y');
    Curl_y = -Bh_xy.*(Dh_x*field_z);
end;

function int_cdu = advect(velx, vely, field);
    int_cdu = tensor2(JtBm_y,JtBm_x,
        tensor2(Jy,Jx,velx).*tensor2(Jy,JD_x,field) +
        tensor2(Jy,Jx,vely).*tensor2(JD_y,Jx,field));
end;

% Solution buffer
P = zeros(nx+1, ny+1);
U = {P, P, P}; U{1}(:,:) = u_ini(x, y);
V = {P, P, P}; V{1}(:,:) = v_ini(x, y);
T = {P, P, P}; T{1}(:,:) = t_ini(x, y);
Vort = P;

% Error buffer
u_err = zeros(1, n_steps);
v_err = zeros(1, n_steps);
t_err = zeros(1, n_steps);
uv_err = zeros(1, n_steps);

% Temporaries
Uh = P; Vh = P; Th = P;
Ut = P; Vt = P;
GradP_x = P; GradP_y = P;

% Advection - Source buffer
AMSx = {P, P, P}; AMSy = {P, P, P}; AMSt = {P, P, P}; 
AMSx{1}(:,:) = (-Bh_xy.*u_src(0,x,y) + advect(U{1},V{1},U{1}))*ht;
AMSy{1}(:,:) = (-Bh_xy.*v_src(0,x,y) + advect(U{1},V{1},V{1}))*ht;
AMSt{1}(:,:) = (-Bh_xy.*t_src(0,x,y) + advect(U{1},V{1},T{1}))*ht;

% Curl of vorticity buffer
CurlVort_x = {P, P, P}; CurlVort_y = {P, P, P};
Vort = compute_vorticity(U{1}, V{1});
[CurlVort_x{1}, CurlVort_y{1}] = curl(ht/Re*Vort);

%% ==================== Setup ends ====================

% Time stepping loop
tic();
for k = 1:n_steps;
    fprintf('Progress: %5.2f%% \r', k / n_steps * 100);
    kc = mod(k, 3) + 1;
    km1 = mod(k - 1, 3) + 1;
    km2 = mod(k - 2, 3) + 1;
    km3 = kc;

    % Set BDFk/EXTk coefficients for first few steps
    if k > 2; 
        a1=ext3(1); a2=ext3(2); a3=ext3(3); 
        b0=bdf3(1); b1=bdf3(2); b2=bdf3(3); b3=bdf3(4);
    elseif k == 2;
        a1=ext2(1); a2=ext2(2); a3=0; 
        b0=bdf2(1); b1=bdf2(2); b2=bdf2(3); b3=0;
    else;
        a1=ext1(1); a2=0; a3=0; 
        b0=bdf1(1); b1=bdf1(2); b2=0; b3=0;
    end;

    % ========== Solution steps start here ===========
    % Step 1: Compute Uh, Vh
    Uh(:,:) = -(
        Bh_xy .* (b1*U{km1} + b2*U{km2} + b3*U{km3})
        + a1*AMSx{km1} + a2*AMSx{km2} + a3*AMSx{km3});
    Vh(:,:) = -(
        Bh_xy .* (b1*V{km1} + b2*V{km2} + b3*V{km3})
        + a1*AMSy{km1} + a2*AMSy{km2} + a3*AMSy{km3});
    Th(:,:) = -(
        Bh_xy .* (b1*T{km1} + b2*T{km2} + b3*T{km3})
        + a1*AMSt{km1} + a2*AMSt{km2} + a3*AMSt{km3});

    % Step 2: Solve for pressure
    Ut(:,:) = Uh - (a1*CurlVort_x{km1} + a2*CurlVort_x{km2} + a3*CurlVort_x{km3});
    Vt(:,:) = Vh - (a1*CurlVort_y{km1} + a2*CurlVort_y{km2} + a3*CurlVort_y{km3});
    P(:,:) = (1/ht)*tensor2(
        Sp_y, Sp_x, Lp_inv.*tensor2(Sp_y', Sp_x', Dh_x'*Ut + Vt*Dh_y));

    GradP_x(:,:) = Dh_x * P; GradP_y(:,:) = P * Dh_y';
    Uh(:,:) = Uh - ht*Bh_xy.*GradP_x;
    Vh(:,:) = Vh - ht*Bh_xy.*GradP_y;

    % Step 3: Viscous solve
    U{kc}(:,:) = tensor2(Su_y, Su_x, Lu_inv.*tensor2(Su_y', Su_x', Uh));
    V{kc}(:,:) = tensor2(Sv_y, Sv_x, Lv_inv.*tensor2(Sv_y', Sv_x', Vh));
    T{kc}(:,:) = tensor2(St_y, St_x, Lt_inv.*tensor2(St_y', St_x', Th));
    
    % ========== Solution steps end here ===========

    % Update advection - source term
    AMSx{kc}(:,:) = ht*(-Bh_xy.*u_src(k*ht,x,y) + advect(U{kc},V{kc},U{kc}));
    AMSy{kc}(:,:) = ht*(-Bh_xy.*v_src(k*ht,x,y) + advect(U{kc},V{kc},V{kc}));
    AMSt{kc}(:,:) = ht*(-Bh_xy.*t_src(k*ht,x,y) + advect(U{kc},V{kc},T{kc}));

    % Update Vort and CurlVort
    Vort = compute_vorticity(U{kc}, V{kc});
    [CurlVort_x{kc}, CurlVort_y{kc}] = curl(ht/Re*Vort);

    % Time-step dependent update
    if k == 1;
        Lu_inv = 1.0 ./ (ht/Re*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf2(1));
        Lv_inv = 1.0 ./ (ht/Re*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf2(1));
        Lt_inv = 1.0 ./ (ht/Pe*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf2(1));
    elseif k == 2;
        Lu_inv = 1.0 ./ (ht/Re*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf3(1));
        Lv_inv = 1.0 ./ (ht/Re*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf3(1));
        Lt_inv = 1.0 ./ (ht/Pe*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf3(1));
    end;    

    % Exact solution: compute L2-error norm
    if exact_u;
        if k < 3; U{kc}(:,:) = u_ext(k*ht,x,y); end;
        u_err(k) = sum(sum(Bh_xy.*(U{kc}-u_ext(k*ht,x,y)).^2));
        u_err(k) = sqrt((0.25/Jac_x/Jac_y)*u_err(k));
    end;

    if exact_v;
        if k < 3; V{kc}(:,:) = v_ext(k*ht,x,y); end;
        v_err(k) = sum(sum(Bh_xy.*(V{kc}-v_ext(k*ht,x,y)).^2));
        v_err(k) = sqrt((0.25/Jac_x/Jac_y)*v_err(k));
    end;

    if exact_t;
        if k < 3; T{kc}(:,:) = t_ext(k*ht,x,y); end;
        t_err(k) = sum(sum(Bh_xy.*(T{kc}-t_ext(k*ht,x,y)).^2));
        t_err(k) = sqrt((0.25/Jac_x/Jac_y)*t_err(k));
    end;

    if exact_u && exact_v;
        uv_err(k) = sum(sum(Bh_xy.*(
            (U{kc}-u_ext(k*ht,x,y)).^2 + (V{kc}-v_ext(k*ht,x,y)).^2
        )));
        uv_err(k) = sqrt((0.25/Jac_x/Jac_y)*uv_err(k));
    end;

    % Plotting
    if k == 1 || k==n_steps || mod(k, ceil(n_steps / 50)) == 0;
        figure(1, 'Units', 'inches', 'Position', [2 2 5 5])
        % quiver(xu,yu,tensor2(Ju_y, Ju_x, U{kc}), tensor2(Ju_y, Ju_x, V{kc}))
        JU = Ju_x * U{kc} * Ju_y';
        mesh(xu,yu,u_ext(k*ht,xu,yu))
        % pbaspect([1 1 1])
        % xlim([bp(1) bp(2)]); ylim([bp(3) bp(4)]);
        % zlim([0 1])
        title(uv_err(k))
        % pause(0.01)
        pause(0.01)
    end;
    
end; % Time loop
et = toc();

U = U{kc}; V = V{kc}; T = T{kc};
soln = {U,V,P,T};
err = {u_err,v_err,t_err,uv_err};
end;
