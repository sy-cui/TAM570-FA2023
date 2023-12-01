addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

N       = [16,16,16];
Tend    = 20; 
Re      = 100; 
Ra      = 10000;
Pr      = 1;
CFL     = 0.1;

nu      = Pr; 
alpha   = 1;

% PCG params
tol=1e-5; max_iter=100; sdim=20;

% Basic setup
nx=N(1);ny=N(2);nz=N(3);
[dim,x,w,D,wm,Jm,Jf] = set_mono_param(N);
[r,s,t] = x{:};
% x = 0.5*(r+1);
% y = 0.5*(s+1); 
% z = 0.5*(t+1); 
x=r;y=s;z=t;
[X,Y,Z] = ndgrid(x,y,z); 
Xf=t3w(Jf,X); Yf=t3w(Jf,Y); Zf=t3w(Jf,Z); 
% [X,Y,Z] = morph_sphere(X,Y,Z,1.0,1.5);
[Rx,G,B,Jac] = geom_elem_3D(X,Y,Z,D,w);
% vol = 4/3*pi*(1.5^3-1)
% scatter3(X(:,:,1),Y(:,:,1),Z(:,:,1)); axis equal
[JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);

% Time steps
dx = min([x(2)-x(1) y(2)-y(1) z(2)-z(1)]);
[dt,nsteps] = cfl_dt(CFL,dx,Tend);

% Operators
Ru=set_restriction(N,{'p','p','d','d','p','p'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D);
Rv=set_restriction(N,{'p','p','d','d','p','p'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D);
Rw=set_restriction(N,{'p','p','d','d','p','p'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D);
Rt=set_restriction(N,{'p','p','d','d','p','p'}); [At,St,Lt]=neumann_op(dim,Rt,w,D);
Rp=set_restriction(N,{'p','p','n','n','p','p'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D);

% Ru=set_restriction(N,{'p','p','p','p','p','p'}); [Au,Su,Lu]=neumann_op(dim,Ru,w,D);
% Rv=set_restriction(N,{'p','p','p','p','p','p'}); [Av,Sv,Lv]=neumann_op(dim,Rv,w,D);
% Rw=set_restriction(N,{'p','p','p','p','p','p'}); [Aw,Sw,Lw]=neumann_op(dim,Rw,w,D);
% Rt=set_restriction(N,{'p','p','p','p','p','p'}); [At,St,Lt]=neumann_op(dim,Rt,w,D);
% Rp=set_restriction(N,{'p','p','p','p','p','p'}); [Ap,Sp,Lp]=neumann_op(dim,Rp,w,D);

% Solution fields
P = zeros(N+1);
U = {P,P,P};
V = {P,P,P};
W = {P,P,P};
T = {P,P,P};

Uhb=P; Vhb=P; Whb=P; Thb=P;  % Inhomogeneous Dirichlet 
Thb(:,1,:) = 1; T{1} = Thb;
% wn = 2*pi; utrans = 1.0; vtrans=1.0;
% U{1} =  sin(wn*X).*cos(wn*Y)+utrans;
% V{1} = -cos(wn*X).*sin(wn*Y)+vtrans; 
epsilon = 1e-3;
% U{1} = epsilon*cos(2*pi*X).*sin(2*pi*Y);
% V{1} = -epsilon*sin(2*pi*X).*sin(2*pi*Y);

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
    Fx=B.*fx(X,Y,Z,km1_time);
    Fy=B.*Ra*Pr.*T{km1};
    Fz=B.*fz(X,Y,Z,km1_time);
    Qt=B.*qt(X,Y,Z,km1_time);

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

    if k < 3;
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
    Uh = Uh;% - viscous_op(Ub,{1,1,1},D,G,B,b0,ndt);
    Vh = Vh;% - viscous_op(Vb,{1,1,1},D,G,B,b0,ndt);
    Wh = Wh;% - viscous_op(Wb,{1,1,1},D,G,B,b0,ndt);
    Th = Th - viscous_op(Tb,{1,1,1},D,G,B,b0,adt);

    % Pressure solve
    rhs = pp_rhs(Ut,Vt,Wt,Rx,D,dt);
    
    % P = tensor3(Sp{3},Sp{2},Sp{1},Lp_inv.*tensor3(Sp{3}',Sp{2}',Sp{1}',rhs));
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
    % mesh(Xf(:,:,2), Yf(:,:,2), t3w(Jf,T{kc})(:,:,2));
    umax = max(max(max(abs(U{kc}))));
    vmax = max(max(max(abs(V{kc}))));
    wmax = max(max(max(abs(W{kc}))));

    % ua = sin(wn*(X-utrans*k*dt)).*cos(wn*(Y-vtrans*k*dt))*exp(-2*ndt*k*wn*wn)+utrans;
    % va =-cos(wn*(X-utrans*k*dt)).*sin(wn*(Y-vtrans*k*dt))*exp(-2*ndt*k*wn*wn)+vtrans;
    % if k == 1 || k == 2;
    %     U{kc} = ua;
    %     V{kc} = va;
    % end;
    % pa =0.25*( cos(wn*2*(X-utrans*k*dt))+cos(wn*2*(Y-vtrans*k*dt)) )*exp(-4*nu*k*dt*wn^2);
    % ue = max(max(max(abs(U{kc} - ua))))
    % ve = max(max(max(abs(V{kc} - va))))
    % disp([it_p,it_u,it_v,it_w,it_t])
    % disp([res_p,res_u,res_v,res_w,res_t])
    
    % if max([umax vmax wmax]) > 5;
    %     error("Unstable");
    % end;

    %% Diagnostics
    if mod(k,20) == 0; 
        figure(1)
        mesh(Xf(:,:,2), Yf(:,:,2), t3w(Jf,T{kc})(:,:,2));
        xlabel("x"); ylabel("y");
        % plot(Jf{1}*x,t3w(Jf,V{kc})(:,ny,nz))
        % hold on
        % plot(Jf{2}*y,t3w(Jf,U{kc})(nx,:,nz))
        % hold off
        % xlim([0 1]); ylim([-1 1])
        % title(sprintf("vmax: %d, umin: %d", max(t3w(Jf,V{kc})(:,ny,nz)), min(t3w(Jf,U{kc})(nx,:,nz))))
        drawnow
    end;


end;

