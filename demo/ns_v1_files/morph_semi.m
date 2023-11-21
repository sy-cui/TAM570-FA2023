function [Xm,Ym] = morph_semi(X,Y);

hdr;

x0=min(min(min(X))); x1=max(max(max(X)));  Lx=x1-x0
y0=min(min(min(Y))); y1=max(max(max(Y)));  Ly=y1-y0

X=(X-x0)/Lx; X=2*(X-.5);  %% On [-1,1]
Y=(Y-y0)/Ly; Y=2*(Y-.5);

R_inner = 1.0;
R_outer = 1.5;
dR = R_outer-R_inner;

Theta = pi*X;
R     = R_inner + dR*0.5*(1+Y);

Xm=R.*sin(Theta);
Ym=R.*cos(Theta);
