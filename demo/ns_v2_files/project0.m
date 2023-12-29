function [Xb,Fb]=project0(Fl,Xk,h0,h1,Bl,Grr,Grs,Gss,Dh);

  n1 = size(Fl,1); n2 = size(Fl,2); n3 = size(Fl,3);
  Fb = reshape(Fl,n1*n2*n3,1);

  alpha = Xk'*Fb;
  Xb    = reshape(Xk*alpha,n1,n2,n3);
  Fb    = Fl-axl(Xb,h0,h1,Bl,Grr,Grs,Gss,Dh);

