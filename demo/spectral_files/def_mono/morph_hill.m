function [Xm,Ym] = morph_hill(X,Y);

hdr;

x0=min(min(X)); x1=max(max(X));  Lx=x1-x0;
y0=min(min(Y)); y1=max(max(Y));  Ly=y1-y0;

X=(X-x0)/Lx; X=2*(X-.5);  %% On [-1,1]
Y=(Y-y0)/Ly; Y=2*(Y-.5);

Xm=X;

Amp = 0.9;
d2 = 0.4^2;
Ym = Y + Amp*0.5*(1-Y).*exp(-X.*X/d2);

Xm = x0 + 0.5*Lx*(Xm+1);
Ym = y0 + 0.5*Ly*(Ym+1);

mesh(Xm,Ym,0*Ym); view(2); axis equal;

