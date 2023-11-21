function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]...
             = set_sem_all(N,Nelx_f,Nelx_b,Nely);

if mod(Nely,2) ~= 1;
   error ("Nely must be odd");
end;

hdr;

Nelx = Nelx_f+Nelx_b+1;
E = Nelx*Nely; Esq = (Nely-1)/2*Nelx + Nelx_f + 1;

% x: Inflow-outflow
% y: Periodic
BC_all = [ 'D' 'N' 'P' 'P' ;     %% U
           'D' 'N' 'P' 'P' ;     %% V
           'N' 'D' 'P' 'P' ;     %% P
           'D' 'N' 'P' 'P' ];    %% T

%% Base Geometry
x0=-0.5-Nelx_f; x1=0.5+Nelx_b; Lx=x1-x0;
y0=-0.5-(Nely-1)/2;    y1=-y0; Ly=y1-y0;

zc = zwuni(Nelx); xc = x0 + Lx*(zc+1)/2;
zc = zwuni(Nely); yc = y0 + Ly*(zc+1)/2;

%% Problem parameters, as function of N

N1=N+1;

[z,w]=zwgll(N);                % Set basic operators
Dh=deriv_mat(z);

[R,S]=ndgrid(z,z);             % Build SEM mesh
X=zeros(N1,E,N1); Y=X;
e=0; 
for ey=1:Nely; for ex=1:Nelx; e=e+1;
    xe0=xc(ex); xe1=xc(ex+1); xeL=xe1-xe0;
    ye0=yc(ey); ye1=yc(ey+1); yeL=ye1-ye0;
    X(:,e,:) = xe0 + xeL*(R+1)/2;
    Y(:,e,:) = ye0 + yeL*(S+1)/2;
end; end;

% [X,Y]=morph_circ(X,Y);         % Morph mesh

[Grr,Grs,Gss,Bl,Xr,Rx,Jac]=geom_elem(X,Y,Dh,w); % Terms for "A"
vol = sum(sum(sum(Bl)))
[Q,glo_num]=set_tp_semq(Nelx,Nely,N);


%% Set masks based on BC_all
   [Mu,Q,glo_num]=set_mask(BC_all(1,:),Nelx,Nely,Q,glo_num); 
   [Mv,Q,glo_num]=set_mask(BC_all(2,:),Nelx,Nely,Q,glo_num); 
   [Mp,Q,glo_num]=set_mask(BC_all(3,:),Nelx,Nely,Q,glo_num); ifnull=1;
   [Mt,Q,glo_num]=set_mask(BC_all(4,:),Nelx,Nely,Q,glo_num); 
   Mu(:,Esq,:)=0; Mu = qqt_op(Mu,'*',glo_num);
   Mv(:,Esq,:)=0; Mv = qqt_op(Mv,'*',glo_num);
   Mp(:,Esq,:)=0; 
   Mt(:,Esq,:)=0; Mt = qqt_op(Mt,'*',glo_num); 
     

%% Set unit-normals for Pressure inhomogeneous BC
   [unxa_v,unya_v] = set_unxy(Mu,X,Y,Xr);

%% Initial conditions
   [U,V,T]=set_ic_bc(0*X,0*X,0*X,Mu,Mv,Mt,X,Y,1);    % 1 --> set IC

%% Diagonal of A matrix (set only after Q is defined)
   dA=diag_sem(Grr,Grs,Gss,Dh); dA=qqt(Q,dA); dA=1./dA;
