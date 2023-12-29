addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

N       = [8,8,8];
Tend    = 100; 
Ra      = 100;
Pr      = 1;
CFL     = 0.01;
l       = 1;
ratio   = 0.5;

nu      = Pr; 
alpha   = 1;
temp_fac = Ra*Pr;

% nu = sqrt(Pr / Ra);
% alpha = 1 / sqrt(Ra * Pr);
% temp_fac = 1;

% PCG params
tol=1e-5; max_iter=1000; sdim=20;

% Basic setup
nx=N(1);ny=N(2);nz=N(3);
[dim,x,w,D,wm,Jm,Jf] = set_mono_param(N,[1,0,0]);
[r,s,t] = x{:};
x=r;y=s;z=t;
[X,Y,Z] = ndgrid(x,y,z); 
[X,Y,Z,Rad,Pol,Azi] = morph_sphere(X,Y,Z,ratio,1.0);
Xf=t3w(Jf,X); Yf=t3w(Jf,Y); Zf=t3w(Jf,Z); 
[Rx,G,B,Jac] = geom_elem_3D(X,Y,Z,D,w);
vol = 4/3*pi*(1-ratio^3)
[JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);

% Time steps
dl = {2*pi,pi,1-ratio};
dx = min([r(2)-r(1) s(2)-s(1) t(2)-t(1)]);
[dt,nsteps] = cfl_dt(CFL,dx,Tend);

% Operators
Ru=set_restriction(N,{'p','p','n','n','d','d'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D,dl);
Rv=set_restriction(N,{'p','p','n','n','d','d'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D,dl);
Rw=set_restriction(N,{'p','p','n','n','d','d'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D,dl);
Rt=set_restriction(N,{'p','p','n','n','d','d'}); [At,St,Lt]=neumann_op(dim,Rt,w,D,dl); 
Rp=set_restriction(N,{'p','p','n','n','n','n'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D,dl);

% Solution fields
P = zeros(N+1);
U = {P,P,P};
V = {P,P,P};
W = {P,P,P};
T = {P,P,P};

Uhb=P; Vhb=P; Whb=P; Thb=P;  % Inhomogeneous Dirichlet 
Thb(:,:,1) = 1; T{1} = Thb;
epsilon = 1e-3;
% Legendre polynomials of cos(\theta)
cP = cos(Pol);
if l==1; Yl0 = sqrt(3/pi)/2*cP; end;
if l==2; Yl0 = sqrt(5/pi)/4*(3*cP.^2-1); end;
if l==3; Yl0 = sqrt(7/pi)/4*(5*cP.^3-3*cP); end;
if l==4; Yl0 = sqrt(9/pi)/16*(35*cP.^4-30*cP.^2+3); end;
if l==5; Yl0 = sqrt(11/pi)/16*(63*cP.^5-70*cP.^3+15*cP); end;
if ~exist("Yl0"); error("Invalid l value"); end;
% U{1} = epsilon.*sin(2*pi*Rad).*Yl0;
% V{1} = epsilon.*sin(2*pi*Rad).*Yl0;
% W{1} = epsilon.*sin(2*pi*Rad).*Yl0;

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

pause

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
    % [P,it_p,res_p,pr] = pcg(t3w(Rp,rhs),Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,max_iter,tol);
    [P,Pk,it_p,res_p] = pcg_prj(rhs,Pk,Pinv_p,Sp,Lp_inv,Rp,D,G,B,0,1,k,sdim,max_iter,tol);
    [gpx,gpy,gpz] = grad_3d(P,Rx,D);
    Uh=Uh-dt*B.*gpx; Vh=Vh-dt*B.*gpy; Wh=Wh-dt*B.*gpz;

    % Viscous solve
    [U{kc},it_u,res_u,~] = pcg(t3w(Ru,Uh),Pinv_u,Su,Lu_inv,Ru,D,G,B,b0,ndt,max_iter,tol);
    [V{kc},it_v,res_v,~] = pcg(t3w(Rv,Vh),Pinv_v,Sv,Lv_inv,Rv,D,G,B,b0,ndt,max_iter,tol);
    [W{kc},it_w,res_w,~] = pcg(t3w(Rw,Wh),Pinv_w,Sw,Lw_inv,Rw,D,G,B,b0,ndt,max_iter,tol);
    [T{kc},it_t,res_t,~] = pcg(t3w(Rt,Th),Pinv_t,St,Lt_inv,Rt,D,G,B,b0,adt,max_iter,tol);

    % disp([it_p it_u it_v it_w it_t])
    U{kc}=U{kc}+Ub; V{kc}=V{kc}+Vb; W{kc}=W{kc}+Wb; T{kc}=T{kc}+Tb;
    umax = max(max(max(abs(U{kc}))));
    vmax = max(max(max(abs(V{kc}))));
    wmax = max(max(max(abs(W{kc}))));
    tmax = max(max(max(abs(T{kc}))));

    disp([it_p,it_u,it_v,it_w,it_t])
    disp([res_p,res_u,res_v,res_w,res_t])
    
    if max(isnan(U{kc})) ~= 0;
        error("Unstable");
    end;

    %% Diagnostics
    if true || mod(k,10) == 0; 
        figure(1)
        % xp = reshape(Xf(:,:,:),[],1); yp = reshape(Yf(:,:,:),[],1); zp = reshape(Zf(:,:,:),[],1);
        % up = reshape(t3w(Jf,T{kc})(:,:,:),[],1);
        % scatter3(xp,yp,zp,[],up); axis equal
        Rf = squeeze(Rad(1,end:-1:1,:));
        Pf = squeeze(Pol(1,end:-1:1,:));
        Tf = squeeze(T{kc}(1,end:-1:1,:));
        contourf(Pf,Rf,Tf)
        % plot(Rf(1,:),Tf(1,:), '-k'); 
        % hold on;
        % plot(Rf(1,:),1./Rf(1,:)-1, '-r');
        % hold off;
        % contourf(squeeze(Xf(:,ny,:)), squeeze(Zf(:,ny,:)), squeeze(t3w(Jf,T{kc})(:,ny,:)));  
        % xlabel("x"); ylabel("z");
        title(num2str([umax, vmax, wmax, tmax]))
        drawnow
    end;


end;

