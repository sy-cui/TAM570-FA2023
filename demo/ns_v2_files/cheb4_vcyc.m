function [Umg,itc,res]=cheb4_vcyc(Fl,nvcyc,k,lam_max,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
    Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,Ue,X,Y);

timec=0;

rl=Fl; Umg=0*Fl; vol = sum(sum(sum(Blc)));
for iter=1:nvcyc % max_iter;

  %% Smoothing step
  [Us]=cheby4(rl,lam_max,0,k,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,Ue,X,Y);
  Umg = Umg+Us;

  %% Coarse-grid correction
  rl = rl-axl(Us,h0,h1,Bl,Grr,Grs,Gss,Dh);
  rc = tensor3(Jc',1,Jc',rl);
  zc = (Blc.*Mc).*(dAc.*rc).^2; zcn=sqrt(sum(sum(sum(zc)))/vol); 
  tolc=0.5*zcn;
  [Ec,itc,res,lam_crs,timec]=...
      pcg_lambda(rc,tolc,10,h0,h1,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,ifnull,timec);
  El = tensor3(Jc,1,Jc,Ec);
  Umg = Umg + El;
  rl  = rl-axl(El,h0,h1,Bl,Grr,Grs,Gss,Dh);
  zl  = (Bl.*M).*(dA.*rl).^2; zn=sqrt(sum(sum(sum(zl)))/vol); 

end;

res = zn;



