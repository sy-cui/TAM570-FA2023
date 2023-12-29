function Xk=project1(Xk,m_max,dX,X,h0,h1,Bl,Grr,Grs,Gss,Dh,M,Q,ifnull);

  N1=size(X,1); E=size(X,2); nl=N1*E*N1;

  m=size(Xk,2);

  if m==m_max; 
     Xk=reshape(X,nl,1)/a_norm(X,h0,h1,Bl,Grr,Grs,Gss,Dh);
  else;
     dW=axl(dX,h0,h1,Bl,Grr,Grs,Gss,Dh);
     dXnrm2 = sum(sum(sum(dW.*dX)));
     dW=reshape(dW,nl,1);
     alpha = Xk'*dW;
     dX=dX-reshape(Xk*alpha,N1,E,N1);
%    anrm = a_norm(dX,h0,h1,Bl,Grr,Grs,Gss,Dh); % Can avoid this extra A*x
     alph2 = sum(alpha.^2);
     anrm  = sqrt(dXnrm2-alph2);
     if anrm>0;
        dX=dX/anrm;
        Xk=[Xk reshape(dX,nl,1)];
     end;
  end;
