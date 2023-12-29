function [x,iter,res]=mg_ksp(Fl,lam_max,tol,max_iter,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                                                   Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,Ue,X,Y);

vol = sum(sum(sum(Blc)));

N1=size(Fl,1); E =size(Fl,2); nl=E*N1*N1;

rho1=1; rtz1=1; r=Fl; x=0*r; p=x;
for iter=1:max_iter;

   [z,itc,rz]=cheb4_vcyc(r,1,5,lam_max,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                               Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,Ue,X,Y);
   rtz0 = rtz1; rtz1=sum(sum(sum(z.*r))); beta=rtz1/rtz0;
   res=sqrt(rtz1);
   if iter > 1 && res < tol; break; end;
   if res == 0; break; end;
   p=z;
   w=axl(p,h0,h1,Bl,Grr,Grs,Gss,Dh);
   pap = sum(sum(sum(p.*w)));
   if iter > 1;
      beta = P'*reshape(w,nl,1);
      p=p-reshape(P*beta,N1,E,N1);
      w=w-reshape(W*beta,N1,E,N1);
      pap=sum(sum(sum(p.*w))); s=1./sqrt(pap);
      P=[P s*reshape(p,nl,1)]; 
      W=[W s*reshape(w,nl,1)];
   else;
      pap=sum(sum(sum(p.*w))); s=1./sqrt(pap);
      P=[s*reshape(p,nl,1)]; 
      W=[s*reshape(w,nl,1)];
   end;
   alpha = s*sum(sum(sum(p.*r)));
   x=x+alpha*p;
   r=r-alpha*w;
end;
