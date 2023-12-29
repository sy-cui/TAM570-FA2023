function [x,iter,res]=mg_ksp(Fl,lam_max,tol,max_iter,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                                                   Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,Ue,X,Y);

format shorte
rl=Fl; z=0*Fl; vol = sum(sum(sum(Blc)));

N1=size(Fl,1); E =size(Fl,2); nl=E*N1*N1;


rho1=1; rtz1=1; r=Fl; x=0*r; p=x;

d=zeros(max_iter,1); l=d; u=d; %% For eigenvalue estimates

hold off;
for iter=1:max_iter;

   [z,itc,rz]=cheb4_vcyc(r,1,6,lam_max,h0,h1,M,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull,...
                               Jc,Mc,Qc,Blc,Grrc,Grsc,Gssc,Dc,dAc,Ue,X,Y);

   rtz0 = rtz1; rtz1=sum(sum(sum(z.*r))); beta=rtz1/rtz0;
   res=sqrt(rtz1);
   if iter > 1 && res < tol; break; end;
   if res == 0; break; end;
   p=z+beta*p;
   if iter==1; Pk=reshape(p,nl,1)/a_norm(p,h0,h1,Bl,Grr,Grs,Gss,Dh); 
        else;  Pk=project1(Pk,40,p,x,h0,h1,Bl,Grr,Grs,Gss,Dh,ifnull); end;
   [dx,r]=project0(r,Pk,h0,h1,Bl,Grr,Grs,Gss,Dh);
   x = x + dx;

end;

if iter<max_iter; iter=iter-1; end;

