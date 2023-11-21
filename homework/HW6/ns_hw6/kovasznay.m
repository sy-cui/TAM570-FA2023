function [U,V] = kovasznay(X,Y,Re);

%  Kovasznay flow, Re=40

lam = Re/2 - sqrt(Re*Re/4 + 4*pi*pi);
U = 1-exp(lam*X).*cos(2*pi*Y);
V = (.5*lam/pi).*exp(lam*X).*sin(2*pi*Y);
