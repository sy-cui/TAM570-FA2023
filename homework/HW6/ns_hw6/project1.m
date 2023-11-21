function Xk=project1(Xk,dX,X,h0,h1,Bl,Grr,Grs,Gss,Dh,M,Q,ifnull);

  N1=size(X,1); E=size(X,2); nl=N1*E*N1;

  m=size(Xk,2);

  m_max=20;

  if m==m_max; 
     Xk=reshape(X,nl,1)/a_norm(X,h0,h1,Bl,Grr,Grs,Gss,Dh);
  else;
     dW=axl(dX,h0,h1,Bl,Grr,Grs,Gss,Dh);
     dW=reshape(dW,nl,1);
     dX=dX-reshape(Xk*(Xk'*dW),N1,E,N1);
     anrm = a_norm(dX,h0,h1,Bl,Grr,Grs,Gss,Dh); % Can avoid this extra A*x
     if anrm>0;
        dX=dX/anrm;
        Xk=[Xk reshape(dX,nl,1)];
     end;
  end;
