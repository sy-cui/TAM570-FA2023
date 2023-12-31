% 2-D SEM convection on periodic domain

%
%    This code supports different Nx and Ny values.
%
%    It also gives the option to have Fourier bases in the 
%    x-direction (ifourier=1), or standard nodal Legendre 
%    bases (ifourier=0).
%
%    As a base case, it solves plane Poiseuille flow.
%
%    It can be modified to solve the Orr-Sommerfeld model 
%    problem.
%


hdr;
x0=0;  x1=2*pi; Lx=x1-x0;   % Domain coordinates and size
y0=-1; y1=1;    Ly=y1-y0;

tfinal = 10;               % Final time

Pr = 10; % Prandtl number
Ra = 3000;
eps = 1.e-4;

nu = Pr;
alpha = 1;

Nx=32; Ny=32; CFL = 0.1;

[Ahx,Bhx,Chx,Dhx,zx,wx] =  semhat(Nx); nhx = Nx+1; Ihx = speye(nhx); 
[Ahy,Bhy,Chy,Dhy,zy,wy] =  semhat(Ny); nhy = Ny+1; Ihy = speye(nhy); 

Mx=ceil(1.2*Nx); [zmx,wmx]=zwgl (Mx); Jx=interp_mat(zmx,zx); Dtx=Jx*Dhx;
My=ceil(1.2*Ny); [zmy,wmy]=zwgl (My); Jy=interp_mat(zmy,zy); Dty=Jy*Dhy;

Nfx = 5+ceil(1.7*Nx); [zfx,w]=zwuni(Nfx); Jfx=interp_mat(zfx,zx); % plotting
Nfy = 5+ceil(1.7*Ny); [zfy,w]=zwuni(Nfy); Jfy=interp_mat(zfy,zy); % plotting

ifourier = 1; if ifourier > 0; % Use nodal Fourier in x (1=yes, 0=no)
 [Ahx,Bhx,Chx,Dhx,zx,wx]= fsemhat(Nx); nhx = Nx+1; Ihx = speye(nhx); 
 Mx=ceil(1.5*Nx);      [zmx,wmx]=zwuni(Mx);  Jx =f_interp_mat(zmx,zx); Dtx=Jx*Dhx;
 Nfx = 5+ceil(1.7*Nx); [zfx,w]=zwuni(Nfx);   Jfx=f_interp_mat(zfx,zx); % plotting
end;

x = x0 + (x1-x0)*(zx+1)/2;  y = y0 + (y1-y0)*(zy+1)/2;  [X,Y]=ndgrid(x,y);  
xf= x0 + (x1-x0)*(zfx+1)/2; yf= y0 + (y1-y0)*(zfy+1)/2; [Xf,Yf]=ndgrid(xf,yf);  

BMxy = (Lx*Ly/4)*wmx*wmy';  %% In shape of domain

Lx2=Lx/2; Lxi=1./Lx2;
Ly2=Ly/2; Lyi=1./Ly2;

Xm=min(min(X)); Ym=min(min(Y));     % For plotting
XM=max(max(X)); YM=max(max(Y));     % For plotting

O=0*X; F=O; G=O; H=O;
U1=O;U2=O;U3=O; V1=O;V2=O;V3=O;
F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
T1=O;T2=O;T3=O; H1=O;H2=O;H3=O;

delta = 0.10*(Lx/2);
Xt = X-.0*(Lx/2);   %  Xt = sin(pi*Xt);
Yt = Y-.0*(Ly/2);   %  Xt = sin(pi*Xt);
R2 = Yt.*Yt + Xt.*Xt;

T  = sin(2*pi*X/Lx).*(1-Y.*Y);
T  = exp(-(R2/(delta*delta)).^1);
%T  = 0.25*(X/pi).*(1-Y);
% T = 00.*X;
% bottomHalf = 1:columns(T)/2;
% T(:,bottomHalf) = 1;
Tb = 00.*Y; Tb(:,1) = 1.0*X(:,1); % Boundary Condition of Temperature

figure(1)
surf(X,Y,T);

Bhh = (wx*wy');
Bxy = (Lx*Ly/4)*Bhh;


RXx=Ihx(2:end-1,:);      RXy=Ihy(2:end-1,:);
RYx=Ihx(2:end-1,:);      RYy=Ihy(2:end-1,:);
RTx=RXx;                 RTy=RYy;
RPx=Ihx;                 RPy=Ihy;  %% Pure Neumann for Pressure

% Set periodicity (here, in x-direction only).

RXx=r_periodic(Nx);  %   RXy=r_periodic(Ny);
RYx=r_periodic(Nx);  %   RYy=r_periodic(Ny);
RTx=r_periodic(Nx);  %   RTy=r_periodic(Ny);
RPx=r_periodic(Nx);  %   RPy=r_periodic(Ny);

Qx = RXx'; Qy = RXy';      % For projection into continuous
Wt = ones(Nx+1,Ny+1);      % space; for plotting purposes only.
Wt = Qx*(Qx'*Wt*Qy)*Qy';
Wt = 1./Wt;


AXx=Lxi*RXx*Ahx*RXx'; BXx=Lx2*RXx*Bhx*RXx'; [SXx,DXx]=gen_eig_ortho(AXx,BXx);
AYx=Lxi*RYx*Ahx*RYx'; BYx=Lx2*RYx*Bhx*RYx'; [SYx,DYx]=gen_eig_ortho(AYx,BYx);
ATx=Lxi*RTx*Ahx*RTx'; BTx=Lx2*RTx*Bhx*RTx'; [STx,DTx]=gen_eig_ortho(ATx,BTx);
APx=Lxi*RPx*Ahx*RPx'; BPx=Lx2*RPx*Bhx*RPx'; [SPx,DPx]=gen_eig_ortho(APx,BPx);

AXy=Lyi*RXy*Ahy*RXy'; BXy=Ly2*RXy*Bhy*RXy'; [SXy,DXy]=gen_eig_ortho(AXy,BXy);
AYy=Lyi*RYy*Ahy*RYy'; BYy=Ly2*RYy*Bhy*RYy'; [SYy,DYy]=gen_eig_ortho(AYy,BYy);
ATy=Lyi*RTy*Ahy*RTy'; BTy=Ly2*RTy*Bhy*RTy'; [STy,DTy]=gen_eig_ortho(ATy,BTy);
APy=Lyi*RPy*Ahy*RPy'; BPy=Ly2*RPy*Bhy*RPy'; [SPy,DPy]=gen_eig_ortho(APy,BPy);

[SoTx,DoTx]=gen_eig_ortho(Lxi*Ahx,Lx2*Bhx); 
[SoTy,DoTy]=gen_eig_ortho(Lyi*Ahy,Ly2*Bhy);

SXx=RXx'*SXx; SXy=RXy'*SXy;  % Prolongate now to save on later multiplies
SYx=RYx'*SYx; SYy=RYy'*SYy;
STx=RTx'*STx; STy=RTy'*STy;
SPx=RPx'*SPx; SPy=RPy'*SPy;

ex=1+0*DXx; ey=1+0*DXy; XLam=ex*DXy' + DXx*ey';
ex=1+0*DYx; ey=1+0*DYy; YLam=ex*DYy' + DYx*ey';
ex=1+0*DTx; ey=1+0*DTy; TLam=ex*DTy' + DTx*ey';
ex=1+0*DPx; ey=1+0*DPy; PLam=ex*DPy' + DPx*ey';

PLam(1,1) = 1;  DPi=1./PLam; DPi(1,1) = 0;  %% Neumann operator for pressure

dx   = min(diff(x)); dy   = min(diff(y)); dx=min(dx,dy);
c    = 1;
if c>0;  dt = CFL*(dx/c);   end; 
if c==0; dt = CFL*dx;       end; 
nsteps = floor(tfinal/dt);
iostep = floor(nsteps/90);
nsteps = iostep*ceil(nsteps/iostep);
dt     = tfinal/nsteps;

tstart=tic; k=0; kk=0; clear tt uu vv ee tnt unt;


U0 = 0*Y + eps * (sin((2*pi*X)/Lx) * cos((pi*Y)/2) + sin((4*pi*X)/Lx) * cos((3*pi*Y)/2));                          %% Set initial conditions
V0 = 0*X + eps * sin((2*pi*X)/Lx);                             %% Set initial conditions

U = U0;
V = V0;

% nsteps = 1;
figure('Position', [100, 100, 1200, 400]); % Adjust the last two values for width and height
for iloop=1:1;
  for istep =1:nsteps; k=k+1; time = k*dt;

     ndt = nu*dt;
     adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     if k<4; DXi=1./(ndt*XLam+b0);DYi=1./(ndt*YLam+b0);DTi=1./(adt*TLam+b0);end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

     Ur = Lxi*BMxy.*(Jx*U*Jy'); Vr = Lyi*BMxy.*(Jx*V*Jy'); 

%%   Set body force, volumetric heating

     if k==1; Q  =  .00*Y; end;  %% No forcing, no heating
     if k==1; FX =  .00*Y;  end;
     if k==1; FY =  Ra * Pr * (T) + .00*Y; end;

%%   Evaluate curl-curl term (to be extrapolated)

     Omega = Lxi*(Dhx*V) - Lyi*(U*Dhy');
     curlcurlX =  Bxy.*(Lyi*(Omega*Dhy'));
     curlcurlY = -Bxy.*(Lxi*(Dhx*Omega));

%%   Compute u-hat and u-tilde

     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advect(U,Ur,Vr,Jx,Jy,Dtx,Dty)+Bxy.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bxy.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);

%%   Compute v-hat and v-tilde

     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advect(V,Ur,Vr,Jx,Jy,Dtx,Dty)+Bxy.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bxy.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);

%%   Compute t-hat and v-tilde

     T3=T2;T2=T1;T1=T;
        H3=H2;H2=H1;H1=-advect(T,Ur,Vr,Jx,Jy,Dtx,Dty)+Bxy.*Q;
        temp = -adt * ((Lxi*Ahx)*Tb*(Ly2*Bhy)' + (Lx2*Bhx)*Tb*(Lyi*Ahy)');
        Th=Bxy.*(b1*T1+b2*T2+b3*T3 - b0*Tb) + (d1*H1+d2*H2+d3*H3) + temp;

%    Pressure correction

     divUT = Lxi*Dhx'*(Ut) + Lyi*(Vt)*Dhy;
     P = (1./dt)*( SPx*(DPi.*(SPx'*divUT*SPy))*SPy');
     dPdx = Lxi*(Dhx*P); dPdy = Lyi*(P*Dhy');

     Uh = Uh - dt*Bxy.*dPdx;
     Vh = Vh - dt*Bxy.*dPdy;

%    Viscous solve

     U = SXx*(DXi.*(SXx'*Uh*SXy))*SXy';  %% Already prolongated
     V = SYx*(DYi.*(SYx'*Vh*SYy))*SYy';  %% to boundary
     T = STx*(DTi.*(STx'*Th*STy))*STy' + Tb;  %%

%    Diagonostics

     umax=max(max(abs(U)));
     vmax=max(max(abs(V)));
     tmax=max(max(abs(T)));

     unt(k)= umax;
     tnt(k)= time;

     if mod(k,10) ==0 || k==1;  kk=kk+1;

       Ut = U; Vt=V;

       vv(kk) = vmax; uu(kk) = umax; tt(kk) = time;

%      Compute vorticity at current step for plotting

       Omega = Lxi*(Dhx*Vt) - Lyi*(Ut*Dhy'); 

       clf;

       Uf=Jfx*Ut*Jfy'; Vf=Jfx*Vt*Jfy'; Tf=Jfx*(Omega)*Jfy';
       %hold off; contour (Xf,Yf,Tf);  axis equal;
       %hold off; surf(X,Y,T);
     %   hold off;  quiver  (Xf,Yf,Uf,Vf,'k-'); axis equal;
       %axis([Xm XM Ym YM ]); axis off;  axis equal;

       subplot(1, 2, 1); quiver  (Xf,Yf,Uf,Vf,'k-'); axis equal; title('UV');

       subplot(1, 2, 2); surf(X,Y,T); title('T');
       xlim([0 2*pi]);
       ylim([-1 1]);
       zlim([0 1]);

       s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(umax) ',   ' ...
                       num2str(vmax) ',  ' num2str(tmax)'.'];
       axes( 'visible', 'off', 'title', s, fs, 12);%(s,fs,12);

       pause(1.0);
     end;
   end;
end;

elapsed_time = toc(tstart);