addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

N       = [32, 16, 16];
Tend    = 5; 
Re      = 100;
Pr      = 1;
CFL     = 0.05;
w2f     = true;
save_dir= 'data';

nu      = 1 / Re; 
alpha   = 1 / (Re * Pr);

% PCG params
tol=1e-5; max_iter=100; sdim=20;

% Basic setup
nx=N(1);ny=N(2);nz=N(3);
[dim,x,w,D,wm,Jm,Jf] = set_mono_param(N);
[r,s,t] = x{:};
x=3*(r+1);y=s;z=t;
[X,Y,Z] = ndgrid(x,y,z); 
[X,Y,Z] = morph_gaussian(X,Y,Z,0.5,0.25,1.5,0.0);
Xf=t3w(Jf,X); Yf=t3w(Jf,Y); Zf=t3w(Jf,Z); 
[Rx,G,B,Jac] = geom_elem_3D(X,Y,Z,D,w);
[JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);

figure(1, 'Units', 'inches', 'Position', [2 2 6 6]); box on;
    hold on
    mesh(Xf(:,:,1), Yf(:,:,1), Zf(:,:,1))
    mesh(Xf(:,:,12), Yf(:,:,12), Zf(:,:,12))
    mesh(Xf(:,:,24), Yf(:,:,24), Zf(:,:,24))
    mesh(Xf(:,:,end), Yf(:,:,end), Zf(:,:,end))
    hold off
    axis equal; view(-45, 15)
    xlabel("$x$", intp, ltx); xlim([0.0 6.0]);
    ylabel("$y$", intp, ltx); ylim([-1.0 1.0]);
    zlabel("$z$", intp, ltx); zlim([-1.0 1.0]);
    set(gca, fn, 'serif', fs, 16, lw, 1.5)
    savefig_png("figures/fph_setup")
pause

if w2f;
    dlmwrite([save_dir, "/X.csv"], reshape(Xf, [], 1), ",");
    dlmwrite([save_dir, "/Y.csv"], reshape(Yf, [], 1), ",");
    dlmwrite([save_dir, "/Z.csv"], reshape(Zf, [], 1), ",");
end;

% Time steps
dl = {x(end)-x(1), y(end)-y(1), z(end)-z(1)};
dx = min([x(2)-x(1) y(2)-y(1) z(2)-z(1)]);
[dt,nsteps] = cfl_dt(CFL,dx,Tend);
nsteps

% Operators
% Inflow-outflow in x; Periodic in y; Rigid wall in z
Ru=set_restriction(N,{'d','n','p','p','d','d'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D,dl);
Rv=set_restriction(N,{'d','n','p','p','d','d'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D,dl);
Rw=set_restriction(N,{'d','n','p','p','d','d'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D,dl);
Rt=set_restriction(N,{'d','n','p','p','d','d'}); [At,St,Lt]=neumann_op(dim,Rt,w,D,dl);
Rp=set_restriction(N,{'n','d','p','p','n','n'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D,dl); ifnull=false;

% Solution fields
P = zeros(N+1);
U = {P,P,P};
V = {P,P,P};
W = {P,P,P};
T = {P,P,P};

Uhb=P; Vhb=P; Whb=P; Thb=P; Phb=P;  % Inhomogeneous Dirichlet 
Uhb(1,:,:) = 1.0; U{1} = Uhb;
% Phb(1,:,:) = 100; P = Phb;

% epsilon = 1e-3;
% U{1} = epsilon*sin(2*pi*X/dl{1}).*cos(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});
% V{1} = -epsilon*cos(2*pi*X/dl{1}).*sin(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});
% W{1} = epsilon*sin(2*pi*X/dl{1}).*sin(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});

% Vorticity and curl of vorticity
Vort_x = P; Vort_y = P; Vort_z = P;
CurlVort_x = {P,P,P};
CurlVort_y = {P,P,P};
CurlVort_z = {P,P,P};

% Advection minus source
ams_x = {P,P,P};
ams_y = {P,P,P};
ams_z = {P,P,P};
ams_t = {P,P,P};

% BDFk/EXTk coefficients
bdf1 = [1 -1];
bdf2 = [3 -4 1] / 2;
bdf3 = [11 -18 9 -2] / 6;
ext1 = [1];
ext2 = [2 -1];
ext3 = [3 -3 1];

% Body force / source term function f(x,t)
fx = @(X,Y,Z,t) 0*X;
fy = @(X,Y,Z,t) 0*X;
fz = @(X,Y,Z,t) 0*X;
qt = @(X,Y,Z,t) 0*X;

% Store A-conjugate basis for pressure
Pk = [];
dA = diag_a_3d(D,G);
Bx = w{2} .* w{3}';

kk = 0;

% Time loop
for k = 1:nsteps;
    fprintf('Progress: %5.2f%% \r', k / nsteps * 100);
    kc  = mod(k,  3)+1;
    km1 = mod(k-1,3)+1;
    km2 = mod(k-2,3)+1;
    km3 = kc;

    % Compute curlcurl and ams from previous step
    [Vort_x,Vort_y,Vort_z]  = curl_3d(U{km1},V{km1},W{km1},Rx,D);
    [CurlVort_x{km1},CurlVort_y{km1},CurlVort_z{km1}]...
                            = curl_3d(Vort_x,Vort_y,Vort_z,Rx,D);
    CurlVort_x{km1} = B.*CurlVort_x{km1};
    CurlVort_y{km1} = B.*CurlVort_y{km1};
    CurlVort_z{km1} = B.*CurlVort_z{km1};

    km1_time = (k-1)*dt;
    Fx=B.*0;
    Fy=B.*0;
    Fz=B.*0;
    Qt=B.*0;

    [Cr,Cs,Ct] = compute_advect_field_3d(U{km1},V{km1},W{km1},Jm,JRx);
    ams_x{km1} = dt*(advect_3d(U{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fx);
    ams_y{km1} = dt*(advect_3d(V{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fy);
    ams_z{km1} = dt*(advect_3d(W{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fz);
    % ams_t{km1} = dt*(advect_3d(T{km1},Cr,Cs,Ct,Jm,JD,Bm) - Qt);

    % Solution steps
    if k == 1;
        a1=ext1(1); a2=0; a3=0; 
        b0=bdf1(1); b1=bdf1(2); b2=0; b3=0;
    elseif k == 2;
        a1=ext2(1); a2=ext2(2); a3=0; 
        b0=bdf2(1); b1=bdf2(2); b2=bdf2(3); b3=0;
    elseif k == 3; 
        a1=ext3(1); a2=ext3(2); a3=ext3(3); 
        b0=bdf3(1); b1=bdf3(2); b2=bdf3(3); b3=bdf3(4);
    end;
    ndt = nu * dt; adt = alpha * dt;

    if k < 4;
        Lu_inv = 1.0 ./ (ndt*Lu + b0);
        Lv_inv = 1.0 ./ (ndt*Lv + b0);
        Lw_inv = 1.0 ./ (ndt*Lw + b0);
        Lt_inv = 1.0 ./ (adt*Lt + b0);
        Lp_inv = 1.0 ./ Lp; if ifnull; Lp_inv(1,1,1)=0; end;
        Pinv_u = 1./t3w(Ru,b0*B+ndt*dA);
        Pinv_v = 1./t3w(Rv,b0*B+ndt*dA);
        Pinv_w = 1./t3w(Rw,b0*B+ndt*dA);
        Pinv_t = 1./t3w(Rt,b0*B+adt*dA);
        Pinv_p = 1./t3w(Rp,dA);
    end;
    
    Ub=U{km1}; Vb=V{km1}; Wb=W{km1}; %Tb=T{km1};
    
    % \hat{u} = -sum(beta_j * u_j + alpha_j * dt * (adv - src)_j)
    Uh = -B.*(b1*U{km1}+b2*U{km2}+b3*U{km3})-(a1*ams_x{km1}+a2*ams_x{km2}+a3*ams_x{km3});
    Vh = -B.*(b1*V{km1}+b2*V{km2}+b3*V{km3})-(a1*ams_y{km1}+a2*ams_y{km2}+a3*ams_y{km3});
    Wh = -B.*(b1*W{km1}+b2*W{km2}+b3*W{km3})-(a1*ams_z{km1}+a2*ams_z{km2}+a3*ams_z{km3});
    % Th = -B.*(b1*T{km1}+b2*T{km2}+b3*T{km3})-(a1*ams_t{km1}+a2*ams_t{km2}+a3*ams_t{km3});

    % \tilde{u} = \hat{u} - nu * dt * curl(vorticity)
    Ut = Uh-ndt*(a1*CurlVort_x{km1}+a2*CurlVort_x{km2}+a3*CurlVort_x{km3});
    Vt = Vh-ndt*(a1*CurlVort_y{km1}+a2*CurlVort_y{km2}+a3*CurlVort_y{km3});
    Wt = Wh-ndt*(a1*CurlVort_z{km1}+a2*CurlVort_z{km2}+a3*CurlVort_z{km3});

    % Include Dirichlet inhomogeneity
    Uh = Uh - viscous_op(Ub,{1,1,1},D,G,B,b0,ndt);
    Vh = Vh - viscous_op(Vb,{1,1,1},D,G,B,b0,ndt);
    Wh = Wh - viscous_op(Wb,{1,1,1},D,G,B,b0,ndt);
    % Th = Th - viscous_op(Tb,{1,1,1},D,G,B,b0,adt);

    % Pressure solve
    flux = X*0; 
    flux(1,:,:) = -squeeze(U{km1}(1,:,:)).*Bx; 
    flux(end,:,:) = squeeze(U{km1}(end,:,:)).*Bx; 
    rhs = pp_rhs(Ut,Vt,Wt,Rx,D,dt) - viscous_op(Phb,{1,1,1},D,G,B,0,1) - (b0/dt)*flux;
    [P,it_p,res_p,pr] = pcg(t3w(Rp,rhs),Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,max_iter,tol);
    % P = t3w(Rp,t3w(Sp,Lp_inv.*t3w(Sp,t3w(Rp,rhs),1)),1);
    P = P + Phb;
    % [P,Pk,it_p,res_p] = pcg_prj(rhs,Pk,Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,k,sdim,max_iter,tol);
    % max(max(max(t3w(Rp,rhs) - viscous_op(pr,Rp,D,G,B,0,1))))
    [gpx,gpy,gpz] = grad_3d(P,Rx,D);
    Uh=Uh-dt*B.*gpx; Vh=Vh-dt*B.*gpy; Wh=Wh-dt*B.*gpz;

    % Viscous solve
    [U{kc},it_u,res_u,~] = pcg(t3w(Ru,Uh),Pinv_u,Su,Lu_inv,Ru,D,G,B,b0,ndt,max_iter,tol);
    [V{kc},it_v,res_v,~] = pcg(t3w(Rv,Vh),Pinv_v,Sv,Lv_inv,Rv,D,G,B,b0,ndt,max_iter,tol);
    [W{kc},it_w,res_w,~] = pcg(t3w(Rw,Wh),Pinv_w,Sw,Lw_inv,Rw,D,G,B,b0,ndt,max_iter,tol);
    % [T{kc},it_t,res_t,~] = pcg(t3w(Rt,Th),Pinv_t,St,Lt_inv,Rt,D,G,B,b0,adt,max_iter,tol);
    % U{kc} = t3w(Ru,t3w(Su,Lu_inv.*t3w(Su,t3w(Ru,Uh),1)),1);
    % V{kc} = t3w(Rv,t3w(Sv,Lv_inv.*t3w(Sv,t3w(Rv,Vh),1)),1);
    % W{kc} = t3w(Rw,t3w(Sw,Lw_inv.*t3w(Sw,t3w(Rw,Wh),1)),1);
    % T{kc} = t3w(Rt,t3w(St,Lt_inv.*t3w(St,t3w(Rt,Th),1)),1);

    U{kc}=U{kc}+Ub; V{kc}=V{kc}+Vb; W{kc}=W{kc}+Wb; %T{kc}=T{kc}+Tb;
    umax = max(max(max(abs(U{kc}))));
    vmax = max(max(max(abs(V{kc}))));
    wmax = max(max(max(abs(W{kc}))));
    % tmax = max(max(max(abs(T{kc}))));
    pmax = max(max(max(abs(P))));

    % disp([it_p,it_u,it_v,it_w])
    % disp([res_p,res_u,res_v,res_w])
    
    if max(isnan(U{kc})) ~= 0;
        error("Unstable");
    end;

    %% Diagnostics
    if mod(k,20) == 0 || k == 1; 
        % figure(1)
        % vmag = sqrt(U{kc}.^2+V{kc}.^2+W{kc}.^2);
        % contourf(squeeze(Xf(:,ny,:)),squeeze(Zf(:,ny,:)),squeeze(t3w(Jf,vmag)(:,ny,:))); axis equal
        % % contourf(squeeze(Xf(:,:,2)),squeeze(Yf(:,:,2)),squeeze(t3w(Jf,V{kc})(:,:,2))); axis equal
        % xlabel("x"); ylabel("z");
        % colorbar;
        % title(num2str([umax vmax wmax pmax]))
        % drawnow

        Ua = t3w(Jf, U{kc});
        Va = t3w(Jf, V{kc});
        Wa = t3w(Jf, W{kc});

        time = k * dt;
        fprintf("Writing to file..., t=%d, kk=%d\n", time, kk);
        dlmwrite([save_dir, sprintf("/U_%d_%d.csv", time, kk)], reshape(Ua,[],1), ",");
        dlmwrite([save_dir, sprintf("/V_%d_%d.csv", time, kk)], reshape(Va,[],1), ",");
        dlmwrite([save_dir, sprintf("/W_%d_%d.csv", time, kk)], reshape(Wa,[],1), ",");

        kk = kk + 1;
    end;
end;

if w2f;
    Ua = t3w(Jf, U{kc});
    Va = t3w(Jf, V{kc});
    Wa = t3w(Jf, W{kc});
    time = k * dt;
    fprintf("Writing to file..., t=%d, kk=%d\n", time, kk);
    dlmwrite([save_dir, sprintf("/U_%d_%d.csv", time, kk)], reshape(Ua,[],1), ",");
    dlmwrite([save_dir, sprintf("/V_%d_%d.csv", time, kk)], reshape(Va,[],1), ",");
    dlmwrite([save_dir, sprintf("/W_%d_%d.csv", time, kk)], reshape(Wa,[],1), ",");

end;

