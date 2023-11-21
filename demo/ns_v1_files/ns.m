clear all; close all;
hdr;    % 2-D SEM multi-element
% close all;

Re=40;
nu=1./Re; alpha=nu;

N=7; 

dt=.025; nsteps=1600;

%% System-solve parameters

ifnull=0; tol=1.e-5; max_iter=140;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]=set_sem_all(N);
[Ue,Ve]=kovasznay(X,Y,Re);

% dx = X(2,1,1) - X(1,1,1);
% CFL = 0.2;
% dt = CFL*dx;

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
Nf = floor(1.2*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);

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

for iloop=1:1;
  for istep =1:nsteps; k=k+1;

%%   Set dt and time
%    [lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5; % Variable dt?
%    dt=ldt_max/lam_max_est;
     time = time+dt;

%%   Set updated BDFk/EXTk coefficients
     ndt = nu*dt; adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=2; a2=-1; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

%%   Set dealiased advecting field
     [Cr,Cs]=set_advect_c(U,V,JM,BMh,Jac,Rx);

%%   Set body force, volumetric heating
     if k==1; QT =  0*Y; end;
     if k==1; FX =  0*Y; end;
     if k==1; FY =  0*Y; end;

%%   Evaluate curl-curl term (to be extrapolated)
     [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Q,Bl,Rx,Dh);


%%   Set Dirichlet conditions onto old fields

     [Ub,Vb,Tb]=set_ic_bc(U,V,T,Mu,Mv,Mt,X,Y,0);         % 0 --> set BC


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
     T3=T2;T2=T1;T1=T;
        H3=H2;H2=H1;H1=-advectl(T,Cr,Cs,JM,DM)+Bl.*QT;
        Th=Bl.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);
        Th=Th-axl(Tb,b0,adt,Bl,Grr,Grs,Gss,Dh);

%    Pressure correction

     divUt = weak_div(Ut,Vt,1.,Rx,Dh)/dt;
%%   Add inhomogeneous Neumann data to divUT, if any. (Eq.(15) in split_slides.pdf)
     b0dt = b0/dt;
     divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub - (1-Mv).*unya_v.*Vb );

%%   Pressure-Poisson solve
     h1=1; h0=0;

%%   Project-out prior solutions, stored in Pk.
     if istep==1;  divUt=divUt-axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh); 
        else;      [P,divUt]=project0(divUt,Pk,h0,h1,Bl,Grr,Grs,Gss,Dh); end;

     [dP,itp,res,lamda_h]=...
         pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,istep);
     P = P+dP;

%%   Save prior solutions, stored in Pk.
     if istep==1; Pk=reshape(P,nl,1)/a_norm(P,h0,h1,Bl,Grr,Grs,Gss,Dh); 
        else;     Pk=project1(Pk,dP,P,h0,h1,Bl,Grr,Grs,Gss,Dh,ifnull); end;

     [dPdx,dPdy]=grad(P,Rx,Dh);
     Uh = Uh - dt*Bl.*dPdx;
     Vh = Vh - dt*Bl.*dPdy;

%    Viscous/diffusive solves (diagonally-preconditioned CG):

%%   Implicit solve - diagonal preconditioner
     dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
     dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+adt*dAU);
     [U,itu,res,lamda_h]=...
        pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
     [V,itv,res,lamda_h]=...
        pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
%    [T,itt,res,lamda_h]=...
%       pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull);
     T=0;

     U=U+Ub;  %% Add back any prescribed Dirichlet conditions
     V=V+Vb;
     T=T+Tb;

     is(istep)=istep; it(istep)=itp;

%    Diagonostics
     if mod(istep,40)==0 || istep<10;  kk=kk+1;

%      disp([itp itu itv itt])

       hold off;

       umax = max(max(max(abs(U))));
       vmax = max(max(max(abs(V))));
       tmax = max(max(max(abs(T))));
       um(kk) = umax; vm(kk) = vmax;
       tm(kk) = tmax; ti(kk) = time;

       Uf = tensor3(Jf,1,Jf,U);  Uf=U;
       Vf = tensor3(Jf,1,Jf,V);  Vf=V;
       Xf = tensor3(Jf,1,Jf,X);  Xf=X;
       Yf = tensor3(Jf,1,Jf,Y);  Yf=Y;
       Tf = tensor3(Jf,1,Jf,T);  Tf=T;

%      if istep>2; Tf=Tf-Tfl; end;
%      Tfl = tensor3(Jf,1,Jf,T);
%      tmax = max(max(max(abs(Tf))));

%      hold off; se_mesh  (Xf,Yf,Tf,s);
%      hold off; se_quiver(Xf,Yf,Uf,Vf,s); axis equal; hold on; 
       drawnow
%      disp([umax vmax tmax])

       dU=Ue-U; dV=Ve-V;
       dU=sqrt(dU.*dU+dV.*dV);
       dinf = max(max(max(dU)));
       ee(kk) = dinf;
%      disp([istep time dinf itp itu itv itt])

       s=['Time: ' num2str(time) ', ' num2str(dinf) ,...
          ', ' int2str([istep itu itv itp])];
       hold off; se_mesh  (Xf,Yf,3*dU./dinf,s);
       ylim([-0.5 0.5]); drawnow
       pause(.04); 
     end;

%    if tmax > 9.9; break; end;
     if umax > 9.9; disp("Unstable!"); break; end;
     if vmax > 9.9; disp("Unstable!"); break; end;

   end;
end;

% save("without_curlcurl.mat", "ti", "ee")
% figure;
% plot(ti,ee,'ko',lw,2);
