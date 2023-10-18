% 2-D NS on periodic domain
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;


alpha = 1.0;
h = 1;
Ly = 2 * h;
Lx = alpha * 2 * pi * h;
tfinal = 60;

U = 1;
Re = 7500;     % Reynolds number
nu = U*h/Re;

CFL = 0.05;

N_arr = [12];
error_buff = 0*N_arr;

for i = 1:length(N_arr);

N = N_arr(i);
[Ah,Bh,Ch,Dh,z,w] = semhat(N);
nh = N+1; Ih = speye(nh); 

M=ceil(1.2*N); [zm,wm]=zwgl(M); Jh=interp_mat(zm,z); Dt=Jh*Dh;
BMxy = (Lx*Ly/4)*wm*wm';

Lx2=Lx/2; Lxi=1./Lx2;
Ly2=Ly/2; Lyi=1./Ly2;

x=Lx*z/2; y=Ly*z/2; [X,Y]=ndgrid(x,y);  

% For plotting:
  Nf=5+ceil(1.7*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z); Xf=Jf*X*Jf'; Yf=Jf*Y*Jf';
  Xm=min(min(X)); Ym=min(min(Y)); XM=max(max(X)); YM=max(max(Y));

% Initialize BDF arrays, etc.
  O=0*X; U=O; V=O; T=O;  F=O; G=O; H=O;
  U1=O;U2=O;U3=O; V1=O;V2=O;V3=O;
  F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
  f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
  T1=O;T2=O;T3=O; H1=O;H2=O;H3=O;

Bhh = (w*w');
Bxy = (Lx*Ly/4)*Bhh;

RXx=Ih(2:end-1,:);       RXy=RXx;
RYx=Ih(2:end-1,:);       RYy=RYx;
RTx=RXx;                 RTy=RYy;
RPx=Ih;                  RPy=RPx;  %% Pure Neumann for Pressure

% Set periodicity (here, in x-direction only).

Nx=N;                    Ny=N; 
RXx=r_periodic(Nx);      %RXy=r_periodic(Ny);
RYx=r_periodic(Nx);      %RYy=r_periodic(Ny);
RTx=r_periodic(Nx);      RTy=r_periodic(Ny);
RPx=r_periodic(Nx);      RPy=r_periodic(Ny);

Qx = RXx'; Qy = RXy';      % For projection into continuous
Wt = ones(Nx+1,Nx+1);      % space; for plotting purposes only.
Wt = Qx*(Qx'*Wt*Qy)*Qy'; Wt = 1./Wt;

AXx = Lxi*RXx*Ah*RXx'; BXx = Lx2*RXx*Bh*RXx'; [SXx,DXx]=gen_eig_decomp(AXx,BXx);
AYx = Lxi*RYx*Ah*RYx'; BYx = Lx2*RYx*Bh*RYx'; [SYx,DYx]=gen_eig_decomp(AYx,BYx);
ATx = Lxi*RTx*Ah*RTx'; BTx = Lx2*RTx*Bh*RTx'; [STx,DTx]=gen_eig_decomp(ATx,BTx);
APx = Lxi*RPx*Ah*RPx'; BPx = Lx2*RPx*Bh*RPx'; [SPx,DPx]=gen_eig_decomp(APx,BPx);
AXy = Lyi*RXy*Ah*RXy'; BXy = Ly2*RXy*Bh*RXy'; [SXy,DXy]=gen_eig_decomp(AXy,BXy);
AYy = Lyi*RYy*Ah*RYy'; BYy = Ly2*RYy*Bh*RYy'; [SYy,DYy]=gen_eig_decomp(AYy,BYy);
ATy = Lyi*RTy*Ah*RTy'; BTy = Ly2*RTy*Bh*RTy'; [STy,DTy]=gen_eig_decomp(ATy,BTy);
APy = Lyi*RPy*Ah*RPy'; BPy = Ly2*RPy*Bh*RPy'; [SPy,DPy]=gen_eig_decomp(APy,BPy);

SXx = RXx'*SXx; SXy = RXy'*SXy;  % Prolongate now to save on later multiplies
SYx = RYx'*SYx; SYy = RYy'*SYy;
STx = RTx'*STx; STy = RTy'*STy;
SPx = RPx'*SPx; SPy = RPy'*SPy;

ex=1+0*DXx; ey=1+0*DXy; XLam=ex*DXy' + DXx*ey';
ex=1+0*DYx; ey=1+0*DYy; YLam=ex*DYy' + DYx*ey';
ex=1+0*DTx; ey=1+0*DTy; TLam=ex*DTy' + DTx*ey';
ex=1+0*DPx; ey=1+0*DPy; PLam=ex*DPy' + DPx*ey';

PLam(1,1) = 1;  DPi=1./PLam; DPi(1,1) = 0;  %% Neumann operator for pressure

%% Estimate CFL and set dt:

   dx   = min(diff(x)); dy   = min(diff(y)); dx=min(dx,dy);
   c = U;
   if c>0;  dt = CFL*(dx/c);   end; 
   if c==0; dt = CFL*dx;       end; 
   nsteps = floor(tfinal/dt);
   iostep = floor(nsteps/90);
   nsteps = iostep*ceil(nsteps/iostep);
   dt     = tfinal/nsteps;

umx=max(max(U))
vmx=max(max(V))
energy_error = zeros(1, nsteps);
% Parabolic profile
   U0 = 1 - Y.^2;
   V0 = 0*Y;
   U = U0 + 1e-5 *sin(2*pi*X/Lx).*cos(2*pi*Y/(2*h));
   V = V0;
   
tstart=tic; k=0; kk=0; clear tt uu tnt unt;
for iloop=1:1;
  for istep =1:nsteps; k=k+1; time = k*dt;

     ndt = nu*dt;
     adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     if k<4; DXi=1./(ndt*XLam+b0); DYi=1./(ndt*YLam+b0); DTi=1./(adt*TLam+b0); end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

     Ur = Lxi*BMxy.*(Jh*U*Jh'); Vr = Lyi*BMxy.*(Jh*V*Jh'); 

     if k==1; Q  =  0.00*Y; end;
     if k==1; FX =  2.0/Re+Y*0; end;
     if k==1; FY =  0.00*Y; end;

     Omega = Lxi*(Dh*V) - Lyi*(U*Dh');
     curlcurlX =  Bxy.*(Lyi*(Omega*Dh'));
     curlcurlY = -Bxy.*(Lxi*(Dh*Omega));

     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advect(U,Ur,Vr,Jh,Dt)+Bxy.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bxy.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);

     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advect(V,Ur,Vr,Jh,Dt)+Bxy.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bxy.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);

     T3=T2;T2=T1;T1=T;H3=H2;H2=H1;H1=-advect(T,Ur,Vr,Jh,Dt)+Bxy.*Q;
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
     Ut = RYx * (U - U0) * RYy';
     Vt = RYx * (V - V0) * RYy';
     energy = sum(sum(tensor2(BYy, BYx, Ut.*Ut + Vt.*Vt)));
     energy_error(istep) = energy;

   %   if k == 1 || mod(k, 1000) == 0;
   %    figure(1)
   %    plot(y, U(ceil(N/2),:))
   %    ylim([-h h])
   %    pause(0.01)
   %    end;

   %   if k == 10000;
   %    break;
   %   end;

   end;
end;
figure(2)
hold on
semilogy([1:nsteps]*dt, energy_error, '-k', lw, 1.5)

l = 0.002234975649;
ae = 1e-9*exp(2*l*[1:nsteps].*dt);
semilogy([1:nsteps]*dt, ae, '-r', lw, 1.5)
hold off


elapsed_time = toc(tstart);

% figure(2)
% quiver(X,Y,U,V)
% xlim([-Lx Lx]/2)
% ylim([-h h])
end;
