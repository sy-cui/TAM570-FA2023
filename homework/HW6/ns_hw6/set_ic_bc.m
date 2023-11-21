function [U,V,T] = set_ic_bc(Uo,Vo,To,Mu,Mv,Mt,X,Y,IC);

%  Kovasznay flow, Re=40
%  Re  = 40;
%  [U,V]=kovasznay(X,Y,Re);

%  Taylor vortices
%  U = -sin(pi*X).*cos(pi*Y);
%  V =  cos(pi*X).*sin(pi*Y);

% Flow past square cylinder
if IC==1;            %% Set Desired IC
   xmax = max(max(max(X))); xmin = min(min(min(X)));
   ymax = max(max(max(Y))); ymin = min(min(min(Y)));
   eps = 1e-1;
   tx = 2*pi/(xmax-xmin)*(X-xmin);
   ty = 2*pi/(ymax-ymin)*(Y-ymin);
   rxy = -(ymax-ymin) / (xmax-xmin);

   U = eps*cos(tx).*sin(ty);
   V = eps*rxy*sin(tx).*cos(ty);
   T = 0 + 0*X;
else;                %% Set Desired BC
   U = 0 + 0*X;
   V = 0 + 0*X;
   T = 0 + 0*X;

   U = Mu.*Uo + (1-Mu).*U;
   V = Mv.*Vo + (1-Mv).*V;
   T = Mt.*To + (1-Mt).*T;
end;
