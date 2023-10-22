addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

%% ==================== Setup ====================
Ra = 1000; Pr = 1; kc = 2*pi;
Tend = 1; CFL = 5e-3;

% Mapping functions from [-1 1]^2
xmap = @(x) (x + 1) * pi / kc; 
ymap = @(y) 0.5*(y + 1);

% Set up basic parameters
ny = 16; nx = ceil(2 * pi / kc * ny);
Jac_x = 0.5 * (xmap(1) - xmap(-1)); % Jacobian in x
Jac_y = 0.5 * (ymap(1) - ymap(-1)); % Jacobian in y

[Ah_x, Bh_x, Ch_x, Dh_x, xi_x, rho_x] = fsemhat(nx);    % Use Fourier in x
[Ah_y, Bh_y, Ch_y, Dh_y, xi_y, rho_y] = semhat(ny);
xh = xmap(xi_x); yh = ymap(xi_y); [x, y] = ndgrid(xh, yh);

% Set up dealiasing
mx = ceil(nx*1.5); my = ceil(ny*1.5);
[eta_x, w_x] = zwgl(mx);
[eta_y, w_y] = zwgl(my);
Jx = f_interp_mat(eta_x, xi_x);
Jy = interp_mat(eta_y, xi_y);

% For plotting
n_uni = 32;
[zu_x, ~] = zwuni(n_uni); zu_x = xmap(zu_x);
[zu_y, ~] = zwuni(n_uni); zu_y = ymap(zu_y);
Ju_x = f_interp_mat(zu_x, xh); Ju_y = interp_mat(zu_y, yh); 
[xu, yu] = ndgrid(zu_x, zu_y);

% Compute time step
hx_min = min(xi_x(2)-xi_x(1)); hy_min = min(xi_y(2)-xi_y(1));
ht = CFL * min(hx_min, hy_min);
n_steps = ceil(Tend / ht);
ht = Tend / n_steps;

% BDFk/EXTk coefficients
bdf1 = [1 -1];
bdf2 = [3 -4 1] / 2;
bdf3 = [11 -18 9 -2] / 6;
ext1 = [1];
ext2 = [2 -1];
ext3 = [3 -3 1];

% Boundary conditions
    % Everything periodic in x
    Ru_x = r_periodic(nx); Rv_x = r_periodic(nx);
    Rp_x = r_periodic(nx); Rt_x = r_periodic(nx);

    Rt_x = speye(nx+1);
    Rp_x = speye(nx+1);
    % Rv_x = speye(nx+1)(2:end-1, :);
    % Ru_x = speye(nx+1)(2:end-1, :);

    % Dirichlet temperature in y
    Rt_y = speye(ny+1)(2:end-1, :);
    % Neumann pressure in y
    Rp_y = speye(ny+1);
    % Homogeneous Dirichlet yvel in y
    Rv_y = speye(ny+1)(2:end-1, :);
    % xvel DD / DN / NN
    Ru_y = speye(ny+1)(2:end-1, :);

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
Bh_xy = (Jac_y*Jac_x)*rho_x.*rho_y';

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
Lu_inv = 1.0 ./ (ht*Pr*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf1(1));
Lv_inv = 1.0 ./ (ht*Pr*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf1(1));
Lt_inv = 1.0 ./ (ht*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf1(1));
Lp_inv = 1.0 ./ (full(diag(Lp_x)) + full(diag(Lp_y))'); Lp_inv(1,1) = 0;

% Helper functions
vorticity = @(velx, vely) Dh_x*vely - velx*Dh_y';
Curl_x = @(field_z) Bh_xy.*(field_z*Dh_y');
Curl_y = @(field_z) -Bh_xy.*(Dh_x*field_z);
advect = @(velx, vely, field) tensor2(JtBm_y,JtBm_x,
    tensor2(Jy,Jx,velx).*tensor2(Jy,JD_x,field) +
    tensor2(Jy,Jx,vely).*tensor2(JD_y,Jx,field));

% Set up inhomogeneous Dirichlet BC for temperature
Tbc_d = zeros(nx+1,ny+1); Tbc_d(:,1) = 1;   % Lower wall T_1 = 1
Tbc_d = sparse(Tbc_d);
H_d = sparse(bdf1(1)*Bh_xy.*Tbc_d+ht*(
    tensor2(Bh_y,Ah_x,Tbc_d)+tensor2(Ah_y,Bh_x,Tbc_d)
));

% Solution buffer
P = zeros(nx+1, ny+1);
U = {P, P, P};
V = {P, P, P};
T = {P, P, P}; T{1}(:,1) = 1;
Vort = P;

% Temporaries
Uh = P; Vh = P; Th = P;
Ut = P; Vt = P;
GradP_x = P; GradP_y = P;

% Advection - Source buffer
AMSx = {P, P, P}; AMSy = {P, P, P}; AMSt = {P, P, P}; 
AMSx{1}(:,:) = advect(U{1},V{1},U{1})*ht;
AMSy{1}(:,:) = (advect(U{1},V{1},V{1}) + Ra*Pr*Bh_xy.*T{1})*ht;
AMSt{1}(:,:) = advect(U{1},V{1},T{1})*ht;

% Curl of vorticity buffer
CurlVort_x = {P, P, P}; CurlVort_y = {P, P, P};
Vort = vorticity(U{1}, V{1});
CurlVort_x{1} = Curl_x(ht*Pr*Vort);
CurlVort_y{1} = Curl_y(ht*Pr*Vort);
CurlVort_x{2} = CurlVort_x{1}; CurlVort_y{2} = CurlVort_y{1}; 
CurlVort_x{3} = CurlVort_x{2}; CurlVort_y{3} = CurlVort_y{2}; 

%% ==================== Setup ends ====================

% Time stepping loop
tic();
for k = 1:n_steps;
    % fprintf('Progress: %5.2f%% \r', k / n_steps * 100);
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
    % Step 1: Compute Uh, Vh, Th
    Uh(:,:) = -(
        Bh_xy .* (b1*U{km1} + b2*U{km2} + b3*U{km3})
        + a1*AMSx{km1} + a2*AMSx{km2} + a3*AMSx{km3});
    Vh(:,:) = -(
        Bh_xy .* (b1*V{km1} + b2*V{km2} + b3*V{km3})
        + a1*AMSy{km1} + a2*AMSy{km2} + a3*AMSy{km3});
    Th(:,:) = -(
        Bh_xy .* (b1*T{km1} + b2*T{km2} + b3*T{km3})
        + a1*AMSt{km1} + a2*AMSt{km2} + a3*AMSt{km3});
    Th(:,:) = Th(:,:) - H_d;

    % Step 2: Solve for pressure
    Ut(:,:) = Uh - ht*Pr*(a1*CurlVort_x{km1} + a2*CurlVort_x{km2} + a3*CurlVort_x{km3});
    Vt(:,:) = Vh - ht*Pr*(a1*CurlVort_y{km1} + a2*CurlVort_y{km2} + a3*CurlVort_y{km3});
    P(:,:) = (1/ht)*tensor2(
        Sp_y, Sp_x, Lp_inv.*tensor2(Sp_y', Sp_x', Dh_x'*Ut + Vt*Dh_y));

    GradP_x(:,:) = Dh_x * P; GradP_y(:,:) = P * Dh_y';
    Uh(:,:) = Uh - ht*Bh_xy.*GradP_x;
    Vh(:,:) = Vh - ht*Bh_xy.*GradP_y;

    % Step 3: Viscous solve
    U{kc}(:,:) = tensor2(Su_y, Su_x, Lu_inv.*tensor2(Su_y', Su_x', Uh));
    V{kc}(:,:) = tensor2(Sv_y, Sv_x, Lv_inv.*tensor2(Sv_y', Sv_x', Vh));
    T{kc}(:,:) = tensor2(St_y, St_x, Lt_inv.*tensor2(St_y', St_x', Th));
    T{kc}(:,:) = T{kc}(:,:) + Tbc_d;

    if sum(sum(isnan(U{kc}))) ~= 0;
        error("Unstable!");
    end;
    
    % ========== Solution steps end here ===========

    % Update advection - source term
    AMSx{kc}(:,:) = ht*advect(U{kc},V{kc},U{kc});
    AMSy{kc}(:,:) = ht*(advect(U{kc},V{kc},V{kc}) + Ra*Pr*Bh_xy.*T{kc});
    AMSt{kc}(:,:) = ht*advect(U{kc},V{kc},T{kc});

    % Update Vort and CurlVort
    Vort(:,:) = vorticity(U{kc}, V{kc});
    CurlVort_x{kc}(:,:) = Curl_x(ht*Pr*Vort);
    CurlVort_y{kc}(:,:) = Curl_y(ht*Pr*Vort);

    % Time-step dependent update
    if k == 1;
        Lu_inv = 1.0 ./ (ht*Pr*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf2(1));
        Lv_inv = 1.0 ./ (ht*Pr*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf2(1));
        Lt_inv = 1.0 ./ (ht*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf2(1));
        H_d = sparse(bdf2(1)*Bh_xy.*Tbc_d+ht*(
            tensor2(Bh_y,Ah_x,Tbc_d)+tensor2(Ah_y,Bh_x,Tbc_d)
        ));
    elseif k == 2;
        Lu_inv = 1.0 ./ (ht*Pr*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf3(1));
        Lv_inv = 1.0 ./ (ht*Pr*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf3(1));
        Lt_inv = 1.0 ./ (ht*(full(diag(Lt_x))+full(diag(Lt_y))') + bdf3(1));
        H_d = sparse(bdf3(1)*Bh_xy.*Tbc_d+ht*(
            tensor2(Bh_y,Ah_x,Tbc_d)+tensor2(Ah_y,Bh_x,Tbc_d)
        ));
    end;    

    %% Plotting
    if k==1 || k==n_steps|| mod(k, ceil(n_steps / 40)) == 0;
        umx = max(max(max(abs(U{kc}))), max(max(abs(V{kc}))))
        figure(1)
        mesh(xu,yu,Ju_x*T{kc}*Ju_y')
        % quiver(xu,yu,tensor2(Ju_y,Ju_x,U{kc}),tensor2(Ju_y,Ju_x,V{kc}))
        axis equal
        view([-120 30])
        pbaspect([1 1 1])
        xlim(xmap([-1 1])); ylim(ymap([-1 1])); zlim([0 1]);
        title(sum(sum(P)))
        pause(0.01)
    end;
    
end; % Time loop
et = toc();
