function [Umg,iter,res]=mg_ksp(Fl,lam_max,tol,max_iter,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                                                   Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc);

format shorte
disp([lam_max tol h0 h1])

rl=Fl; Umg=0*Fl; vol = sum(sum(sum(Blc)));


for iter=1:max_iter;

  %% Smoothing step
  omega=0; if omega > 0;
    [Us]=smoother(rl,lam_max,omega,5,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
  else;
%   [Us]=cheby4(rl,lam_max,omega,6,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
    [Us,iterc,res,lam_crs]=...
      pcg_lambda(rl,tol,500,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
  end;
  Umg = Umg+Us;

  %% Coarse-grid correction
  rl = rl-axl(Us,h0,h1,Bl,Grr,Grs,Gss,Dh);
% rc = tensor3(Jc',1,Jc',rl);
% zc = (Blc.*Mc).*(dAc.*rc).^2; zcn=sqrt(sum(sum(sum(zc)))/vol); 
% tolc=0.1*zcn;
% [Ec,iterc,res,lam_crs]=...
%     pcg_lambda(rc,tolc,100,h0,h1,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,ifnull);
% El = tensor3(Jc,1,Jc,Ec);
% Umg = Umg + El;
% rl  = rl-axl(El,h0,h1,Bl,Grr,Grs,Gss,Dh);

  zl = (Bl.*M).*(dA.*rl).^2; zn=sqrt(sum(sum(sum(zl)))/vol); 

  if zn < tol; break; end;

  format shorte;
  disp([iter zn tol])

end;

res = zn;
iter=min(iter,max_iter);


