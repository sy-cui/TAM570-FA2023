function w = ax(u,b0,nu,Rr,Rs,Bb,Grr,Grs,Gss,Dhr,Dhs,ifnull);

%
%   A = D^T G D
%
%
%   Here, we assume that u is in the restricted space
%

w  = Rr'*u*Rs;
w  = axb(w,b0,nu,Bb,Grr,Grs,Gss,Dhr,Dhs,ifnull);
w  = Rr*w*Rs';

