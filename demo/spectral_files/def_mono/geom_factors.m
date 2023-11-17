function [G,B,Xr,Rx,Jac]=geom_factors(X,Y,Dr,Ds,wr,ws,ifourierx,ifouriery);

hdr;

nr = length(wr);  Nr=nr-1;
ns = length(ws);  Ns=ns-1;

G =zeros(nr,ns,2,2);
Xr=zeros(nr,ns,2,2);
Rx=zeros(nr,ns,2,2);

%  Compute dx_i / dr_j

%
%  First, subtract off mean (nonperiodic component)
%

if ifourierx > 0;                              %% This is NOT perfect...
   [zr,w]=zwuni(Nr);
   X0=sum(X(1,:))/ns; X1=sum(X(end,:))/ns;     %% It's really a kludge
   Y0=sum(Y(1,:))/ns; Y1=sum(Y(end,:))/ns;     %% It's really a kludge
   dXdrb = (X1-X0)/2; X = X-(X0 + (X1-X0)*(1+zr)/2);
   dYdrb = (Y1-Y0)/2; Y = Y-(Y0 + (Y1-Y0)*(1+zr)/2);
else;
   dXdrb = 0; dYdrb = 0;
end;

if ifouriery > 0;
   [zs,w]=zwuni(Ns);
   X0=sum(X(:,1))/nr; X1=sum(X(:,end))/nr;
   Y0=sum(Y(:,1))/nr; Y1=sum(Y(:,end))/nr;
   dXdsb = (X1-X0)/2; X = X-(X0 + (X1-X0)*(1+zs')/2);
   dYdsb = (Y1-Y0)/2; Y = Y-(Y0 + (Y1-Y0)*(1+zs')/2);
else;
   dXdsb = 0; dYdsb = 0;
end;




Xr(:,:,1,1) = Dr*X  + dXdrb;                          % dx/dr: 
Xr(:,:,1,2) = X*Ds' + dXdsb;                          % dx/ds
Xr(:,:,2,1) = Dr*Y  + dYdrb;                          % dy/dr: 
Xr(:,:,2,2) = Y*Ds' + dYdsb;                          % dy/ds

Jac = Xr(:,:,1,1).*Xr(:,:,2,2) - Xr(:,:,1,2).*Xr(:,:,2,1); 
Ji  = 1./Jac;

Rx(:,:,1,1) =  Ji.*Xr(:,:,2,2);                       % dr/dx
Rx(:,:,1,2) = -Ji.*Xr(:,:,1,2);                       % dr/dy
Rx(:,:,2,1) = -Ji.*Xr(:,:,2,1);                       % ds/dx
Rx(:,:,2,2) =  Ji.*Xr(:,:,1,1);                       % ds/dy

B = Jac.*(wr*ws');  % nr x ns diagonal mass matrix on Omega-hat, incl. Jac

G(:,:,1,1)=B.*(Rx(:,:,1,1).*Rx(:,:,1,1)+Rx(:,:,1,2).*Rx(:,:,1,2));% G_rr=sum_j dr/dx_j * dr/dx_j
G(:,:,1,2)=B.*(Rx(:,:,1,1).*Rx(:,:,2,1)+Rx(:,:,1,2).*Rx(:,:,2,2));% G_rs=sum_j dr/dx_j * ds/dx_j
G(:,:,2,1)=G(:,:,1,2);                                            % G_sr=G_rs - symmetric tensor
G(:,:,2,2)=B.*(Rx(:,:,2,1).*Rx(:,:,2,1)+Rx(:,:,2,2).*Rx(:,:,2,2));% G_rr=sum_j dr/dx_j * dr/dx_j

