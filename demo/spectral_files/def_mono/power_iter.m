hdr;    % 2-D SEM convection on periodic domain

%  This version addresses nonperiodic geometry


%% Geometry

imorph=2;

   x0 = -1.5*pi;  x1=-x0; Lx=x1-x0;   % Domain coordinates and size
   y0 =  1.0;      y1=1.5; Ly=y1-y0;

%  x0 = -1;  x1 = 1; Lx=x1-x0;   % Domain coordinates and size
%  y0 = -1;  y1 = 1; Ly=y1-y0;   % Domain coordinates and size


Nx = 127; ifourierx=1; % Nodal Fourier in x (1=yes, 0=no)
Ny = 040; ifouriery=0; % Nodal Fourier in y (1=yes, 0=no)

if ifourierx>0;[zx,wx]=zwuni(Nx);else;[zx,wx]=zwgll(Nx);end;x=x0+Lx*(zx+1)/2;
if ifouriery>0;[zy,wy]=zwuni(Ny);else;[zy,wy]=zwgll(Ny);end;y=y0+Ly*(zy+1)/2;



[Ahx,Bhx,Chx,Dhx,Ihx,Jx,zx,wx,Jfx,zfx]= lfsemhat(Nx,ifourierx);  % Set basic operators
[Ahy,Bhy,Chy,Dhy,Ihy,Jy,zy,wy,Jfy,zfy]= lfsemhat(Ny,ifouriery);  % Set basic operators
xf=x0+Lx*(zfx+1)/2; yf=y0+Ly*(zfy+1)/2;

%%% FOR PLOTTING THE MESH, use fine points _along_ the line, nodes across lines

[X,Y]  =ndgrid(x,y);  [Xf,Yf]=ndgrid(xf,yf);
[X1,Y1]=ndgrid(x,yf); [X2,Y2]=ndgrid(xf,y);

if imorph>0; [X,Y]=morph_hill(X,Y);     [Xf,Yf]=morph_hill(Xf,Yf); 
             [X1,Y1]=morph_hill(X1,Y1); [X2,Y2]=morph_hill(X2,Y2); end;

if imorph>1; [X,Y]=morph_circ(X,Y);     [Xf,Yf]=morph_circ(Xf,Yf);
             [X1,Y1]=morph_circ(X1,Y1); [X2,Y2]=morph_circ(X2,Y2); end;

plot(X2,Y2,'k-',lw,1.3',X1',Y1','k-',lw,1.3); axis equal; axis off;

[G,Bb,Xr,Rx,Jac]=geom_factors(X,Y,Dhx,Dhy,wx,wy,ifourierx,ifouriery);


%% Set boundary restriction operators in "x" and "y" directions (really, "r" and "s" directions)

if ifourierx>0; 
   RXx=r_periodic(Nx); RYx=RXx; RTx=RXx; RPx=RXx; 
else;
   RXx=Ihx(2:end-1,:);RYx=Ihx(2:end-1,:);RTx=RXx;RPx=Ihx;  %% Pure Neumann for Pressure
%  RXx=r_periodic(Nx);RYx=RXx;RTx=RXx;RPx=RXx; 
end;

if ifouriery>0; 
   RXy=r_periodic(Ny); RYy=RXy; RTy=RXy; RPy=RXy; 
else;
   RXy=Ihy(2:end-1,:); RYy=Ihy(2:end-1,:); RTy=RXy; RPy=Ihy;  %% Pure Neumann for Pressure
 % RXy=r_periodic(Ny); RYy=RXy; RTy=RXy; RPy=RXy; 
end;

B=RXx*Bb*RXy';
Bi=1./B;


%% Set 1D (r,s) operators for FDM preconditioner

[SXx,DXx]=set_1d_ops(Ahx,Bhx,RXx,Lx); [SYx,DYx]=set_1d_ops(Ahx,Bhx,RYx,Lx);
[SPx,DPx]=set_1d_ops(Ahx,Bhx,RPx,Lx); [STx,DTx]=set_1d_ops(Ahx,Bhx,RTx,Lx);

[SXy,DXy]=set_1d_ops(Ahy,Bhy,RXy,Ly); [SYy,DYy]=set_1d_ops(Ahy,Bhy,RYy,Ly);
[SPy,DPy]=set_1d_ops(Ahy,Bhy,RPy,Ly); [STy,DTy]=set_1d_ops(Ahy,Bhy,RTy,Ly);


ex=1+0*DXx; ey=1+0*DXy; XLam=ex*DXy' + DXx*ey';
ex=1+0*DYx; ey=1+0*DYy; YLam=ex*DYy' + DYx*ey';
ex=1+0*DTx; ey=1+0*DTy; TLam=ex*DTy' + DTx*ey';
ex=1+0*DPx; ey=1+0*DPy; PLam=ex*DPy' + DPx*ey';

PLam(1,1) = 1;  DPi=1./PLam; DPi(1,1) = 0;  %% Neumann operator for pressure (UNLESS OUTFLOW!)

%
%  Here, we perform our first numerical experiment: Condition number of FDM vs A
%


nr = size(zx,1);
ns = size(zy,1);
nb = nr*ns;       %% Total number of points, including boundaries

Grr = G(:,:,1,1);
Grs = G(:,:,1,2);
Gss = G(:,:,2,2);
Da  = diag_a(Grr,Grs,Gss,Dhx,Dhy);
if ifourierx>0; Da(1,:) = 0.5*Da(2,:); Da(end,:)=Da(1,:); end;
if ifouriery>0; Da(:,1) = 0.5*Da(:,2); Da(:,end)=Da(:,1); end;

F = 1+.1*rand(nr,ns);
F = rand(nr,ns);
Ue= 1+1*rand(nr,ns);
Ue= l2(Ue,RXx,RXy,Bb); % Project into correct space


b0  = 0; nu = 1;

Fb = axb(Ue,b0,nu,Bb,Grr,Grs,Gss,Dhx,Dhy);
Di = nu.*Da + b0; Di=RXx*Di*RXy'; Di=1./Di;

D2=sqrt(Di);

tol = 1.e-9; max_iter = 800; ifnull=0;


format shorte;

x=RXx*Ue*RXy';
xb=RXx'*x*RXy;
for iter=1:max_iter;

    xnorm=sqrt(sum(sum(x)));
    xb=xb/xnorm; x=x/xnorm;
    zb=axb(xb,b0,nu,Bb,Grr,Grs,Gss,Dhx,Dhy);
    z =Di.*(RXx*zb*RXy');
    lam = sum(sum(z));
    x=z; xb=RXx'*x*RXy;
    if mod(iter,80)==0; disp([iter lam]), end;

end;

