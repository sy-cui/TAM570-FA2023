function [x,iter,residual]=cg_h(Fb,tol,max_iter,b0,nu,Rr,Rs,Bb,Grr,Grs,Gss,Dhr,Dhs,ifnull,Di);


% Solve H u = f via CG,  H = b0*B + nu*A

%
%   A = D^T G D
%

% Fb  - rhs, in grid form (full grid)
% Gij - geometric factors, in grid form
% Dhr - derivative in Omega-hat, r direction
% Dhs - derivative in Omega-hat, s direction
% Rr  - 1D restriction matrix, r direction
% Rs  - 1D restriction matrix, s direction
%
% ifnull > 0 --- Implies that H has a null space of 1
%


%  In the code below, use sum(sum(sum(z.*r))) for z'*r.
%  In this way, we can keep z and r in "grid shape" either for 2D or 3D


rho1=1; r=Rr*Fb*Rs'; x=0*r; p=x;

for iter=1:max_iter;

%  z=r;     rho0 = rho1; rho1=sum(sum(sum(z.*r))); beta=rho1/rho0;
   z=Di.*r; rho0 = rho1; rho1=sum(sum(sum(z.*r))); beta=rho1/rho0; % Diagonal preconditioner
   residual=sqrt(rho1);
   if residual < tol; break; end;
   p=z+beta*p;
   w=ax(p,b0,nu,Rr,Rs,Bb,Grr,Grs,Gss,Dhr,Dhs,ifnull);

   gamma = sum(sum(sum(p.*w))); alpha = rho1/gamma;
   x=x+alpha*p;
   r=r-alpha*w;

%  disp([iter residual tol])

end;

if iter < max_iter; iter=iter-1; end;

x = Rr'*x*Rs;  %% Prolongate back to full grid



