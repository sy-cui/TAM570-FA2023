hdr;    % 2-D SEM multi-element
% close all;

Re=40; 
nu=1./Re; alpha=nu;

N=7; 

dt=.0100; nsteps=1500;

%% System-solve parameters

ifnull=0; tol=1.e-5; max_iter=140;
ifnull=0; tol=1.e-2; max_iter=140;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,...
   unxa_v,unya_v,BC_all,dA,Dc,Jc,Grrc,Grsc,Gssc,Blc,Qc,glo_numc,Mc,dAc]=set_sem_mg_all(N);

[Ue,Ve] = kovasznay(X,Y,Re);  %% Exact solution for Kovasznay test case



E =size(U,2);
N1=N+1;
nl=N1*E*N1;

%% Set dealiasing operators, JM,DM,BMh
Nd = floor(1.5*N);
[zd,wd]=zwgl(Nd); Bd=diag(wd);
JM=interp_mat(zd,z);
DM=deriv_mat (zd);
BMh=tensor3(Bd*JM,1,Bd*JM,1+0*Jac); %% effectively, Bd*JM*Jac*JM'*Bd'

%% Set plotting interpolation operator
Nf = floor(2.1*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);

%% Automatic timestep size selection -- FOR LATER

% [lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5;
% lam_max_est=max(1e-5,lam_max_est);
% dt=ldt_max/lam_max_est;

% dt=.01;
% Tfinal = 4*pi; nsteps = ceil(Tfinal/dt)  
% dt = Tfinal/nsteps;

%% Initialize BDFk/EXTk arrays
O=0*X; F=O; G=O; H=O;
U1=O;U2=O;U3=O; V1=O;V2=O;V3=O;
F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
T1=O;T2=O;T3=O; H1=O;H2=O;H3=O;




%%%%% TIME STEPPING LOOP %%%%%%

P=O;  %% Initialize pressure to zero

kk=0; k=0; time=0;

t_tot=0;
t_adv=0;
t_prs=0;
t_vis=0;


for iloop=1:1;
  for istep =1:nsteps; k=k+1; t0=tic;

%%   Set dt and time
%    [lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5; % Variable dt?
%    dt=ldt_max/lam_max_est;
     time = time+dt;

%%   Set updated BDFk/EXTk coefficients
     ndt = nu*dt; adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

%%   Set dealiased advecting field
     [Cr,Cs,t_adv]=set_advect_c(U,V,JM,BMh,Jac,Rx,t_adv);

%%   Set body force, volumetric heating
     if k==1; QT =  0*Y; end;  %% No forcing, no heating
     if k==1; FX =  0*Y; end;
     if k==1; FY =  0*Y; end;

%%   Evaluate Bl*(curl-curl) term (to be extrapolated)
     [curlcurlX,curlcurlY,Omega]=curlcurl_0(U,V,Q,Bl,Rx,Dh);


%%   Set Dirichlet conditions onto old fields

     [Ub,Vb,Tb]=set_ic_bc(U,V,T,Mu,Mv,Mt,X,Y,istep);  % 0 --> set IC


%%   Compute u-hat and u-tilde
     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advectl(U,Cr,Cs,JM,DM)+Bl.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bl.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);
        Uh=Uh-axl(Ub,b0,ndt,Bl,Grr,Grs,Gss,Dh);

%%   Compute v-hat and v-tilde
     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advectl(V,Cr,Cs,JM,DM)+Bl.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bl.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);
        Vh=Vh-axl(Vb,b0,ndt,Bl,Grr,Grs,Gss,Dh);

%%   Compute t-hat
%    T3=T2;T2=T1;T1=T;
%       H3=H2;H2=H1;H1=-advectl(T,Cr,Cs,JM,DM)+Bl.*QT;
%       Th=Bl.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);
%       Th=Th-axl(Tb,b0,adt,Bl,Grr,Grs,Gss,Dh);

%%   Pressure correction

     divUt = weak_div(Ut,Vt,1.,Rx,Dh)/dt;
%%   Add inhomogeneous Neumann data to divUT, if any. (Eq.(15) in split_slides.pdf)
     b0dt = b0/dt;
     divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub - (1-Mv).*unya_v.*Vb );

%%   Pressure-Poisson solve
     h1=1; h0=0;

%%   Project-out prior solutions, stored in Pk.
     if istep==1;  divUt=divUt-axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh); 
        else;      [P,divUt]=project0(divUt,Pk,h0,h1,Bl,Grr,Grs,Gss,Dh); end;

%%   Iterative solve for pressure
     ifnull = istep; % just for diagnostics
     if istep==1; [dP,itp,res,lambda_h,t_prs]=pcg_lambda(divUt,tol,2*max_iter,...
                    h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,t_prs); 
        lambda_h = 1.1*lambda_h
        itc = 0;
     else;
       [dP,itp,itc,res,t_prs]=mg_ksp(divUt,lambda_h,tol,max_iter,...
                    h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                    Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,dP,X,Y,t_prs); end;
     P = P+dP;

%%   Save prior solutions, stored in Pk.
     if istep==1; Pk=reshape(P,nl,1)/a_norm(P,h0,h1,Bl,Grr,Grs,Gss,Dh); 
        else;     Pk=project1(Pk,30,dP,P,h0,h1,Bl,Grr,Grs,Gss,Dh,ifnull); end;

%%   Update velocity field by subtracting pressure correction
     [dPdx,dPdy]=grad(P,Rx,Dh);
     Uh = Uh - dt*Bl.*dPdx;
     Vh = Vh - dt*Bl.*dPdy;

%    Viscous/diffusive solves (diagonally-preconditioned CG):

%%   Implicit solve - diagonal preconditioner
     dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
     dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+adt*dAU);
     [U,itu,res,lamda_h,t_vis]=...
        pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull,t_vis);
     [V,itv,res,lamda_h,t_vis]=...
        pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull,t_vis);
%    [T,itt,res,lamda_h,t_vis]=...
%       pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull,t_vis);
     T=Omega;


     U=U+Ub;  %% Add back any prescribed Dirichlet conditions
     V=V+Vb;
     T=T+Tb;

     is(istep)=istep; it(istep)=itp;

%    Diagonostics

     if istep==1;
        x0=10.1; y0=0.1;
        R=(X-x0).^2+(Y-y0).^2;
        rmin = min(min(min(R)));
        imin = find(R<1.0001*rmin);
     end;
     st_time(istep) = time;
     st_vely(istep) = V(imin);

     t1=toc(t0);
     t_tot = t_tot + t1;

     if mod(istep,40)==0 || istep<3;  kk=kk+1;


       hold off;

       umax = max(max(max(abs(U))));
       vmax = max(max(max(abs(V))));
       tmax = max(max(max(abs(T))));
       um(kk) = umax; vm(kk) = vmax;
       tm(kk) = tmax; ti(kk) = time;

       Xf = tensor3(Jf,1,Jf,X);
       Yf = tensor3(Jf,1,Jf,Y);
       Pf = tensor3(Jf,1,Jf,P);
       Uf = tensor3(Jf,1,Jf,U);
       Vf = tensor3(Jf,1,Jf,V);
       Tp = Pf + .5*(Uf.*Uf+Vf.*Vf);
       Er = sqrt((Ue-U).^2 + (Ve-V).^2);
       Ef = tensor3(Jf,1,Jf,Er); emx=max(max(max(Ef)));

       dU=sqrt(V.*V);
       vinf = max(max(max(dU)));

       s=['Time: ' num2str([time emx  t_tot t_adv t_prs t_vis]) ...
          ', ' int2str([istep itu itv itp itc])];
       hold off; se_quiver (X,Y,U,V,s); 
%      hold off; se_contour(Xf,Yf,Tp,s,20); 
       hold on;  se_mesh(Xf,Yf,3*Ef/emx,s,20); 
                 axis equal; 
       xlabel('X',fs,20); ylabel('Y',fs,20); zlabel('Normalized Error',fs,20);
       drawnow

       pause(.04); 
     end;

%    if tmax > 9.9; break; end;
     if umax > 19.9; break; end;
     if vmax > 19.9; break; end;

   end;
end;

figure;
%   plot(is,it,'ko',lw,2); axis([0 nsteps 0 150]);
semilogy(is,it,'ko',lw,2); axis([0 nsteps 1 150]);
xlabel('Step Number',fs,14);
ylabel('Number of Pressure Iterations per Step',fs,14);
title(['Pressure Iteration Count, Kovasznay, ' int2str([E N])],fs,16);
