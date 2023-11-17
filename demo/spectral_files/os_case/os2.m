% 2-D SEM convection on periodic domain

hdr;

for N=32;

Lx=2*pi;       % Domain size
Ly=2; 

Pr = 5.e-4;   % Prandtl number
Re = 7500,    % Reynolds number
Re = 5000,    % Reynolds number
Re = 6000,    % Reynolds number
Re = 5500,    % Reynolds number
Re = 5800,    % Reynolds number
Re = 7500,    % Reynolds number
nu = 1/Re;
alpha=nu/Pr;


tmax = 1;

tfinal = 100; 
CFL    = 1.0;


[Ah,Bh,Ch,Dh,z,w] = semhat(N);
nh = N+1; Ih = speye(nh); 

M=ceil(1.2*N); [zm,wm]=zwgl(M); Jh=interp_mat(zm,z); Dt=Jh*Dh;
BMxy = (Lx*Ly/4)*wm*wm';



Lx2=Lx/2; Lxi=1./Lx2;
Ly2=Ly/2; Lyi=1./Ly2;

x=Lx*z/2; y=Ly*z/2; [X,Y]=ndgrid(x,y);  

Nf = 5+ceil(1.7*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);   % For plotting
Xf=Jf*X*Jf'; Yf=Jf*Y*Jf';           % For plotting
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

Bhh = (w*w');
Bxy = (Lx*Ly/4)*Bhh;



RXx=Ih(2:end-1,:);       RXy=RXx;
RYx=Ih(2:end-1,:);       RYy=RYx;
RTx=RXx;                 RTy=RYy;
RPx=Ih;                  RPy=RPx;  %% Pure Neumann for Pressure

% Set periodicity (here, in x-direction only).

Nx=N;                %   Ny=N; 
RXx=r_periodic(Nx);  %   RXy=r_periodic(Ny);
RYx=r_periodic(Nx);  %   RYy=r_periodic(Ny);
RTx=r_periodic(Nx);  %   RTy=r_periodic(Ny);
RPx=r_periodic(Nx);  %   RPy=r_periodic(Ny);

Qx = RXx'; Qy = RXy';      % For proejction into continuous
Wt = ones(Nx+1,Nx+1);      % space; for plotting purposes only.
Wt = Qx*(Qx'*Wt*Qy)*Qy';
Wt = 1./Wt;


AXx = Lxi*RXx*Ah*RXx'; BXx = Lx2*RXx*Bh*RXx'; [SXx,DXx]=gen_eig_ortho(AXx,BXx);
AYx = Lxi*RYx*Ah*RYx'; BYx = Lx2*RYx*Bh*RYx'; [SYx,DYx]=gen_eig_ortho(AYx,BYx);
ATx = Lxi*RTx*Ah*RTx'; BTx = Lx2*RTx*Bh*RTx'; [STx,DTx]=gen_eig_ortho(ATx,BTx);
APx = Lxi*RPx*Ah*RPx'; BPx = Lx2*RPx*Bh*RPx'; [SPx,DPx]=gen_eig_ortho(APx,BPx);

AXy = Lyi*RXy*Ah*RXy'; BXy = Ly2*RXy*Bh*RXy'; [SXy,DXy]=gen_eig_ortho(AXy,BXy);
AYy = Lyi*RYy*Ah*RYy'; BYy = Ly2*RYy*Bh*RYy'; [SYy,DYy]=gen_eig_ortho(AYy,BYy);
ATy = Lyi*RTy*Ah*RTy'; BTy = Ly2*RTy*Bh*RTy'; [STy,DTy]=gen_eig_ortho(ATy,BTy);
APy = Lyi*RPy*Ah*RPy'; BPy = Ly2*RPy*Bh*RPy'; [SPy,DPy]=gen_eig_ortho(APy,BPy);


SXx = RXx'*SXx; SXy = RXy'*SXy;  % Prolongate now to save on later multiplies
SYx = RYx'*SYx; SYy = RYy'*SYy;
STx = RTx'*STx; STy = RTy'*STy;
SPx = RPx'*SPx; SPy = RPy'*SPy;

ex=1+0*DXx; ey=1+0*DXy; XLam=ex*DXy' + DXx*ey';
ex=1+0*DYx; ey=1+0*DYy; YLam=ex*DYy' + DYx*ey';
ex=1+0*DTx; ey=1+0*DTy; TLam=ex*DTy' + DTx*ey';
ex=1+0*DPx; ey=1+0*DPy; PLam=ex*DPy' + DPx*ey';

PLam(1,1) = 1;  DPi=1./PLam; DPi(1,1) = 0;  %% Neumann operator for pressure

dx   = min(diff(x)); dy   = min(diff(y)); dx=min(dx,dy);
c    = 2;
if c>0;  dt = CFL*(dx/c);   end; 
if c==0; dt = CFL*dx;       end; 
nsteps = floor(tfinal/dt);
iostep = floor(nsteps/90);
nsteps = iostep*ceil(nsteps/iostep);
dt     = tfinal/nsteps;


U0 = 1-Y.*Y;                          %% Set initial conditions
V0 = 0*X;                             %% Set initial conditions

small = .00005; Xs = 2*pi*X/Lx; Ys = 2*pi*Y/Ly;
U = U0 + small*sin(Ys).*cos(Xs);
V = V0 - 0*small*cos(Ys).*sin(Xs);


umx=max(max(U))
vmx=max(max(V))

tstart=tic; k=0; kk=0; clear tt uu vv ee tnt unt;
for iloop=1:1;
  for istep =1:nsteps; k=k+1; time = k*dt;

     ndt = nu*dt;
     adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     if k<4; DXi=1./(ndt*XLam+b0); DYi=1./(ndt*YLam+b0); DTi=1./(adt*TLam+b0); end;
%    if k == floor(nsteps/2); 
%       Re = 2*Re/3,  % Reynolds number
%       Re = 5770;    % Reynolds number
%       nu = 1/Re;
%       ndt = nu*dt;
%       DXi=1./(ndt*XLam+b0); DYi=1./(ndt*YLam+b0); DTi=1./(adt*TLam+b0);
%       FX =  (2/Re) + 0*Y;
%    end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

     Ur = Lxi*BMxy.*(Jh*U*Jh'); Vr = Lyi*BMxy.*(Jh*V*Jh'); 

     if k==1; Q  =  .00*Y; end;  %% No forcing, no heating
     if k==1; FX =  (2/Re) + 0*Y;  end;
     if k==1; FY =  .00*Y; end;


     Omega = Lxi*(Dh*V) - Lyi*(U*Dh');
     curlcurlX =  Bxy.*(Lyi*(Omega*Dh'));
     curlcurlY = -Bxy.*(Lxi*(Dh*Omega));

     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advects(U,Ur,Vr,Jh,Dt)+Bxy.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bxy.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);

     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advects(V,Ur,Vr,Jh,Dt)+Bxy.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bxy.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);

     T3=T2;T2=T1;T1=T;H3=H2;H2=H1;H1=-advects(T,Ur,Vr,Jh,Dt)+Bxy.*Q;
        Th=Bxy.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);

%    pressure correction
     divUT = Lxi*Dh'*(Ut) + Lyi*(Vt)*Dh;
     P = (1./dt)*( SPx*(DPi.*(SPx'*divUT*SPy))*SPy');
     dPdx = Lxi*(Dh*P); dPdy = Lyi*(P*Dh');

     Uh = Uh - dt*Bxy.*dPdx;
     Vh = Vh - dt*Bxy.*dPdy;

%    viscous solve
     U = SXx*(DXi.*(SXx'*Uh*SXy))*SXy';  %% Already prolongated to boundary
     V = SYx*(DYi.*(SYx'*Vh*SYy))*SYy';
     T = STx*(DTi.*(STx'*Th*STy))*STy';

%    Diagonostics
     umax=max(max(abs(U)));
     vmax=max(max(abs(V)));
     tmax=max(max(abs(T)));

     unt(k)= umax;
     tnt(k)= time;

%    if vmax > 30; break; end;

     if mod(k,100) ==0 || k==1;  kk=kk+1;
       vv(kk) = vmax; uu(kk) = umax; tt(kk) = time;
       s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(umax) ',   ' ...
                       num2str(vmax) ',  ' num2str(istep)'.'];
       Ut = U-U0; Vt=V-V0;
       e2 = sqrt( sum(sum( (Bh*(Ut.*Ut)*Bh + Bh*(Vt.*Vt)*Bh)/4.)) );
       ee(kk)=e2;
       Omega = Lxi*(Dh*Vt) - Lyi*(Ut*Dh');
%      Omega = Wt.*(Qx*(Qx'*Omega*Qy)*Qy');    % Project H^1_0

       Uf=Jf*Ut*Jf'; Vf=Jf*Vt*Jf'; Tf=Jf*(Omega)*Jf';
       hold off; quiver  (X,Y,U,V,'k-'); axis equal; 
       hold off;   contour (Xf,Yf,Tf); axis equal;% axis square;
       hold on;  quiver  (Xf,Yf,Uf,Vf,'k-'); axis equal; 
       axis([Xm XM Ym YM ]); axis off; %axis square;
       title(s,fs,12); drawnow; 
     end;

   end;
end;

elapsed_time = toc(tstart);

figure
plot(tt,vv,'b-',lw,2,tt,ee,'r-',lw,2); axis square;
title('V_{max} vs time',fs,14)
xlabel('time',fs,12); ylabel('V_{max}',fs,12);

Ne(N)=max(unt);
Nn(N)=N;

end;

% tt1=tt;vv1=vv;uu1=uu;ee1=ee;

