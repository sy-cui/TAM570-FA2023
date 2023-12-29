addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

N       = [12,12,12];
Tend    = 200; 
Ra      = 1725;
Pr      = 1;
CFL     = 0.03;
kc      = 3.117/sqrt(2);
sigma   = 1e-4;

nu      = Pr; 
alpha   = 1;
temp_fac = Ra*Pr;

% nu = sqrt(Pr / Ra);
% alpha = 1 / sqrt(Ra * Pr);
% temp_fac = 1;

% PCG params
tol=1e-5; max_iter=10; sdim=20;

% Basic setup
nx=N(1);ny=N(2);nz=N(3);
[dim,x,w,D,wm,Jm,Jf] = set_mono_param(N);
[r,s,t] = x{:};
x = pi/kc*(r+1);
y = pi/kc*(s+1); 
z = 0.5*(t+1);
% x=r;y=s;z=t;
[X,Y,Z] = ndgrid(x,y,z); 
Xf=t3w(Jf,X); Yf=t3w(Jf,Y); Zf=t3w(Jf,Z); 
% [X,Y,Z] = morph_sphere(X,Y,Z,1.0,1.5);
[Rx,G,B,Jac] = geom_elem_3D(X,Y,Z,D,w);
% vol = 4/3*pi*(1.5^3-1)
[JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);

% Time steps
dl = {x(end)-x(1), y(end)-y(1), z(end)-z(1)};
dx = min([x(2)-x(1) y(2)-y(1) z(2)-z(1)]);
[dt,nsteps] = cfl_dt(CFL,dx,Tend);

% Operators
Ru=set_restriction(N,{'p','p','p','p','d','d'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D,dl);
Rv=set_restriction(N,{'p','p','p','p','d','d'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D,dl);
Rw=set_restriction(N,{'p','p','p','p','d','d'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D,dl);
Rt=set_restriction(N,{'p','p','p','p','d','d'}); [At,St,Lt]=neumann_op(dim,Rt,w,D,dl);
Rp=set_restriction(N,{'p','p','p','p','n','n'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D,dl);

% Ru=set_restriction(N,{'p','p','p','p','p','p'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D,dl);
% Rv=set_restriction(N,{'p','p','p','p','p','p'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D,dl);
% Rw=set_restriction(N,{'p','p','p','p','p','p'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D,dl);
% Rt=set_restriction(N,{'p','p','p','p','p','p'}); [At,St,Lt]=neumann_op(dim,Rt,w,D,dl);
% Rp=set_restriction(N,{'p','p','p','p','p','p'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D,dl);

% Solution fields
P = zeros(N+1);
U = {P,P,P};
V = {P,P,P};
W = {P,P,P};
T = {P,P,P};

Uhb=P; Vhb=P; Whb=P; Thb=P;  % Inhomogeneous Dirichlet 
Thb(:,:,1) = 1; T{1} = Thb;
% wn = 2*pi; utrans = 1.0; vtrans=1.0;
% U{1} =  sin(wn*X).*cos(wn*Z)+utrans;
% W{1} = -cos(wn*X).*sin(wn*Z)+vtrans; 
epsilon = 1e-3;
U{1} = epsilon*sin(2*pi*X/dl{1}).*cos(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});
V{1} = -epsilon*cos(2*pi*X/dl{1}).*sin(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});
W{1} = epsilon*sin(2*pi*X/dl{1}).*sin(2*pi*Y/dl{2}).*sin(2*pi*Z/dl{3});

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

KE = zeros(1, nsteps);
growth_rate = sigma + 1;

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
    Fz=B.*temp_fac.*T{km1};
    Qt=B.*0;

    [Cr,Cs,Ct] = compute_advect_field_3d(U{km1},V{km1},W{km1},Jm,JRx);
    ams_x{km1} = dt*(advect_3d(U{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fx);
    ams_y{km1} = dt*(advect_3d(V{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fy);
    ams_z{km1} = dt*(advect_3d(W{km1},Cr,Cs,Ct,Jm,JD,Bm) - Fz);
    ams_t{km1} = dt*(advect_3d(T{km1},Cr,Cs,Ct,Jm,JD,Bm) - Qt);

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
        Lp_inv = 1.0 ./ Lp;  Lp_inv(1,1,1) = 0;
        Pinv_u = 1./t3w(Ru,b0*B+ndt*dA);
        Pinv_v = 1./t3w(Rv,b0*B+ndt*dA);
        Pinv_w = 1./t3w(Rw,b0*B+ndt*dA);
        Pinv_t = 1./t3w(Rt,b0*B+adt*dA);
        Pinv_p = 1./t3w(Rp,dA);
    end;
    
    Ub=U{km1}; Vb=V{km1}; Wb=W{km1}; Tb=T{km1};
    
    % \hat{u} = -sum(beta_j * u_j + alpha_j * dt * (adv - src)_j)
    Uh = -B.*(b1*U{km1}+b2*U{km2}+b3*U{km3})-(a1*ams_x{km1}+a2*ams_x{km2}+a3*ams_x{km3});
    Vh = -B.*(b1*V{km1}+b2*V{km2}+b3*V{km3})-(a1*ams_y{km1}+a2*ams_y{km2}+a3*ams_y{km3});
    Wh = -B.*(b1*W{km1}+b2*W{km2}+b3*W{km3})-(a1*ams_z{km1}+a2*ams_z{km2}+a3*ams_z{km3});
    Th = -B.*(b1*T{km1}+b2*T{km2}+b3*T{km3})-(a1*ams_t{km1}+a2*ams_t{km2}+a3*ams_t{km3});

    % \tilde{u} = \hat{u} - nu * dt * curl(vorticity)
    Ut = Uh-ndt*(a1*CurlVort_x{km1}+a2*CurlVort_x{km2}+a3*CurlVort_x{km3});
    Vt = Vh-ndt*(a1*CurlVort_y{km1}+a2*CurlVort_y{km2}+a3*CurlVort_y{km3});
    Wt = Wh-ndt*(a1*CurlVort_z{km1}+a2*CurlVort_z{km2}+a3*CurlVort_z{km3});

    % Include Dirichlet inhomogeneity
    Uh = Uh - viscous_op(Ub,{1,1,1},D,G,B,b0,ndt);
    Vh = Vh - viscous_op(Vb,{1,1,1},D,G,B,b0,ndt);
    Wh = Wh - viscous_op(Wb,{1,1,1},D,G,B,b0,ndt);
    Th = Th - viscous_op(Tb,{1,1,1},D,G,B,b0,adt);

    % Pressure solve
    rhs = pp_rhs(Ut,Vt,Wt,Rx,D,dt);
    [P,it_p,res_p,pr] = pcg(t3w(Rp,rhs),Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,max_iter,tol);
    % P = t3w(Rp,t3w(Sp,Lp_inv.*t3w(Sp,t3w(Rp,rhs),1)),1);
    % [P,Pk,it_p,res_p] = pcg_prj(rhs,Pk,Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,k,sdim,max_iter,tol);

    % max(max(max(t3w(Rp,rhs) - viscous_op(t3w(Rp,rhs),Rp,D,G,B,0,1))))
    [gpx,gpy,gpz] = grad_3d(P,Rx,D);
    Uh=Uh-dt*B.*gpx; Vh=Vh-dt*B.*gpy; Wh=Wh-dt*B.*gpz;

    % Viscous solve
    [U{kc},it_u,res_u,~] = pcg(t3w(Ru,Uh),Pinv_u,Su,Lu_inv,Ru,D,G,B,b0,ndt,max_iter,tol);
    [V{kc},it_v,res_v,~] = pcg(t3w(Rv,Vh),Pinv_v,Sv,Lv_inv,Rv,D,G,B,b0,ndt,max_iter,tol);
    [W{kc},it_w,res_w,~] = pcg(t3w(Rw,Wh),Pinv_w,Sw,Lw_inv,Rw,D,G,B,b0,ndt,max_iter,tol);
    [T{kc},it_t,res_t,~] = pcg(t3w(Rt,Th),Pinv_t,St,Lt_inv,Rt,D,G,B,b0,adt,max_iter,tol);
    % U{kc} = t3w(Ru,t3w(Su,Lu_inv.*t3w(Su,t3w(Ru,Uh),1)),1);
    % V{kc} = t3w(Rv,t3w(Sv,Lv_inv.*t3w(Sv,t3w(Rv,Vh),1)),1);
    % W{kc} = t3w(Rw,t3w(Sw,Lw_inv.*t3w(Sw,t3w(Rw,Wh),1)),1);
    % T{kc} = t3w(Rt,t3w(St,Lt_inv.*t3w(St,t3w(Rt,Th),1)),1);

    % disp([it_p it_u it_v it_w it_t])
    U{kc}=U{kc}+Ub; V{kc}=V{kc}+Vb; W{kc}=W{kc}+Wb; T{kc}=T{kc}+Tb;
    umax = max(max(max(abs(U{kc}))));
    vmax = max(max(max(abs(V{kc}))));
    wmax = max(max(max(abs(W{kc}))));
    tmax = max(max(max(abs(T{kc}))));

    % ua = sin(wn*(X-utrans*k*dt)).*cos(wn*(Z-vtrans*k*dt))*exp(-2*ndt*k*wn*wn)+utrans;
    % % va =-cos(wn*(X-utrans*k*dt)).*sin(wn*(Y-vtrans*k*dt))*exp(-2*ndt*k*wn*wn)+vtrans;
    % va = 0*X;
    % wa =-cos(wn*(X-utrans*k*dt)).*sin(wn*(Z-vtrans*k*dt))*exp(-2*ndt*k*wn*wn)+vtrans;
    % if k == 1 || k == 2;
    %     U{kc} = ua;
    %     W{kc} = wa;
    % end;
    % % pa =0.25*( cos(wn*2*(X-utrans*k*dt))+cos(wn*2*(Y-vtrans*k*dt)) )*exp(-4*nu*k*dt*wn^2);
    % ue = max(max(max(abs(U{kc} - ua))))
    % ve = max(max(max(abs(V{kc} - va))))
    % we = max(max(max(abs(W{kc} - wa))))
    % disp([it_p,it_u,it_v,it_w,it_t])
    % disp([res_p,res_u,res_v,res_w,res_t])
    
    if max(isnan(U{kc})) ~= 0;
        error("Unstable");
    end;

    KE(k) = sum(sum(sum(0.5 * B.*(U{kc}.^2 + V{kc}.^2 + W{kc}.^2))));
    if k > 1; growth_rate = abs(KE(k) - KE(k-1)) / (dt* KE(k)); end;
    if k*dt > 1 && growth_rate < sigma;
        fprintf("Reached steady state at T = %d\n. KE is %d\n", dt * k, KE(k));
        break;
    end;

    %% Diagnostics
    if mod(k,40) == 0; 
        % figure(1)
        % contourf(squeeze(Xf(:,ny,:)), squeeze(Zf(:,ny,:)), squeeze(t3w(Jf,T{kc})(:,ny,:)));  
        % % plot(Jf{3}*z,squeeze(t3w(Jf,T{kc})(nx,ny,:)))
        % xlabel("x"); ylabel("y");
        % title(num2str([umax, vmax, wmax, tmax growth_rate]))
        % plot(Jf{1}*x,t3w(Jf,V{kc})(:,ny,nz))
        % hold on
        % plot(Jf{2}*y,t3w(Jf,U{kc})(nx,:,nz))
        % hold off
        % xlim([0 1]); ylim([-1 1])
        % title(sprintf("vmax: %d, umin: %d", max(t3w(Jf,V{kc})(:,ny,nz)), min(t3w(Jf,U{kc})(nx,:,nz))))
        drawnow
    end;


end;
KE = KE(1:k);

