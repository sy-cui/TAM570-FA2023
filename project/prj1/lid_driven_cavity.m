addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

%% ==================== Setup ====================
Re = 100; Tend = 30; CFL = 0.2; U0 = 1;
nu = U0 / Re;

% Mapping functions from [-1 1]^2
xmap = @(x) 0.5*(x + 1); 
ymap = @(y) 0.5*(y + 1);

% Set up basic parameters
nx = 64; ny = nx;
Lx = xmap(1) - xmap(-1); Ly = ymap(1) - ymap(-1);
Jac_x = 0.5 * Lx; Jac_y = 0.5 * Ly;

[Ah_x, Bh_x, Ch_x, Dh_x, xi_x, rho_x] = semhat(nx);    % Use Fourier in x
[Ah_y, Bh_y, Ch_y, Dh_y, xi_y, rho_y] = semhat(ny);
xh = xmap(xi_x); yh = ymap(xi_y); [x, y] = ndgrid(xh, yh);

% Set up dealiasing
mx = ceil(nx*1.5); my = ceil(ny*1.5);
[eta_x, w_x] = zwgl(mx);
[eta_y, w_y] = zwgl(my);
Jx = interp_mat(eta_x, xi_x);
Jy = interp_mat(eta_y, xi_y);

% For plotting
n_uni = 64;
[zu_x, ~] = zwuni(n_uni); zu_x = xmap(zu_x);
[zu_y, ~] = zwuni(n_uni); zu_y = ymap(zu_y);
Ju_x = interp_mat(zu_x, xh); Ju_y = interp_mat(zu_y, yh); 
[xu, yu] = ndgrid(zu_x, zu_y);

% Compute time step
hx_min = min(xh(2)-xh(1)); hy_min = min(yh(2)-yh(1));
ht = CFL * min(hx_min, hy_min) / U0;
n_steps = ceil(Tend / ht);
ht = Tend / n_steps;
nu_ht = nu * ht;

% BDFk/EXTk coefficients
bdf1 = [1 -1];
bdf2 = [3 -4 1] / 2;
bdf3 = [11 -18 9 -2] / 6;
ext1 = [1];
ext2 = [2 -1];
ext3 = [3 -3 1];

% Boundary conditions
Ru_x = speye(nx+1)(2:end-1, :); Ru_y = speye(ny+1)(2:end-1, :);
Rv_x = speye(nx+1)(2:end-1, :); Rv_y = speye(ny+1)(2:end-1, :);
Rp_x = speye(nx+1); Rp_y = speye(ny+1);

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
Bh_xy = (Jac_y*Jac_x)*rho_x.*rho_y';

JtBm_x=Jx'*diag(w_x)*Jac_x; JtBm_y=Jy'*diag(w_y)*Jac_y; 
JD_x = Jx*Dh_x; JD_y = Jy*Dh_y;

% Eigen decomposition
[Su_x,Lu_x]=gen_eig_decomp(Au_x,Bu_x); [Su_y,Lu_y]=gen_eig_decomp(Au_y,Bu_y);
[Sv_x,Lv_x]=gen_eig_decomp(Av_x,Bv_x); [Sv_y,Lv_y]=gen_eig_decomp(Av_y,Bv_y);
[Sp_x,Lp_x]=gen_eig_decomp(Ap_x,Bp_x); [Sp_y,Lp_y]=gen_eig_decomp(Ap_y,Bp_y);
Su_x = Ru_x'*Su_x; Su_y = Ru_y'*Su_y;
Sv_x = Rv_x'*Sv_x; Sv_y = Rv_y'*Sv_y;
Sp_x = Rp_x'*Sp_x; Sp_y = Rp_y'*Sp_y;
Lu_inv = 1.0 ./ (nu_ht*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf1(1));
Lv_inv = 1.0 ./ (nu_ht*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf1(1));
Lp_inv = 1.0 ./ (full(diag(Lp_x)) + full(diag(Lp_y))'); Lp_inv(1,1) = 0;

% Helper functions
vorticity = @(velx, vely) Dh_x*vely - velx*Dh_y';
Curl_x = @(field_z) Bh_xy.*(field_z*Dh_y');
Curl_y = @(field_z) -Bh_xy.*(Dh_x*field_z);
advect = @(velx, vely, field) tensor2(JtBm_y,JtBm_x,
    tensor2(Jy,Jx,velx).*tensor2(Jy,JD_x,field) +
    tensor2(Jy,Jx,vely).*tensor2(JD_y,Jx,field));

% Set up inhomogeneous Dirichlet BC for the lid
Ubc_d = zeros(nx+1,ny+1); Ubc_d(:,end) = U0; 
Ubc_d = sparse(Ubc_d);
H_d = sparse(bdf1(1)*Bh_xy.*Ubc_d+nu_ht*(
    tensor2(Bh_y,Ah_x,Ubc_d)+tensor2(Ah_y,Bh_x,Ubc_d)
));

% Solution buffer
P = zeros(nx+1, ny+1);
U = {P, P, P}; U{1} = U{1} + Ubc_d;
V = {P, P, P};
Vort = P;

% Temporaries
Uh = P; Vh = P;
Ut = P; Vt = P;
GradP_x = P; GradP_y = P;

% Advection - Source buffer
AMSx = {P, P, P}; AMSy = {P, P, P};
AMSx{1}(:,:) = advect(U{1},V{1},U{1})*ht;
AMSy{1}(:,:) = advect(U{1},V{1},V{1})*ht;

% Curl of vorticity buffer
CurlVort_x = {P, P, P}; CurlVort_y = {P, P, P};
Vort = vorticity(U{1}, V{1});
CurlVort_x{1} = Curl_x(nu_ht*Vort);
CurlVort_y{1} = Curl_y(nu_ht*Vort);
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
    % Step 1: Compute Uh, Vh
    Uh(:,:) = -(
        Bh_xy .* (b1*U{km1} + b2*U{km2} + b3*U{km3})
        + a1*AMSx{km1} + a2*AMSx{km2} + a3*AMSx{km3});
    Vh(:,:) = -(
        Bh_xy .* (b1*V{km1} + b2*V{km2} + b3*V{km3})
        + a1*AMSy{km1} + a2*AMSy{km2} + a3*AMSy{km3});
    Uh(:,:) = Uh - H_d;


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
    U{kc}(:,:) = U{kc} + Ubc_d;

    if sum(sum(isnan(U{kc}))) ~= 0;
        error("Unstable!");
    end;
    
    % ========== Solution steps end here ===========

    % Update advection - source term
    AMSx{kc}(:,:) = ht*advect(U{kc},V{kc},U{kc});
    AMSy{kc}(:,:) = ht*advect(U{kc},V{kc},V{kc});

    % Update Vort and CurlVort
    Vort(:,:) = vorticity(U{kc}, V{kc});
    CurlVort_x{kc}(:,:) = Curl_x(nu_ht*Vort);
    CurlVort_y{kc}(:,:) = Curl_y(nu_ht*Vort);

    % Time-step dependent update
    if k == 1;
        Lu_inv = 1.0 ./ (nu_ht*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf2(1));
        Lv_inv = 1.0 ./ (nu_ht*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf2(1));
        H_d = sparse(bdf2(1)*Bh_xy.*Ubc_d+nu_ht*(
            tensor2(Bh_y,Ah_x,Ubc_d)+tensor2(Ah_y,Bh_x,Ubc_d)
        ));
    elseif k == 2;
        Lu_inv = 1.0 ./ (nu_ht*(full(diag(Lu_x))+full(diag(Lu_y))') + bdf3(1));
        Lv_inv = 1.0 ./ (nu_ht*(full(diag(Lv_x))+full(diag(Lv_y))') + bdf3(1));
        H_d = sparse(bdf3(1)*Bh_xy.*Ubc_d+nu_ht*(
            tensor2(Bh_y,Ah_x,Ubc_d)+tensor2(Ah_y,Bh_x,Ubc_d)
        ));
    end;    

    %% Plotting
    if k==1 || k==n_steps|| mod(k, ceil(n_steps / 40)) == 0;
        centerV = tensor2(Ju_y,Ju_x,V{kc})(:,n_uni/2);
        [Vmax, imax] = max(centerV);
        [Vmin, imin] = min(centerV);

        figure(1)
        % quiver(xu,yu,tensor2(Ju_y,Ju_x,U{kc}),tensor2(Ju_y,Ju_x,V{kc}))
        plot(xu, tensor2(Ju_y,Ju_x,V{kc})(:,n_uni/2))
        axis equal
        % view([-120 30])
        pbaspect([1 1 1])
        xlim(xmap([-1 1])); ylim([-0.5 0.5])
        title(sprintf("T:%d, Vmax at %d: %d, Vmin at %d: %d", 
            k*ht, xu(imax), Vmax, xu(imin), Vmin)); 
        drawnow
    end;
    
end; % Time loop
et = toc();
