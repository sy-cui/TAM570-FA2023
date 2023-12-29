function [U,V,T] = set_ic_bc(Uo,Vo,To,Mu,Mv,Mt,X,Y,istep);

%  Kovasznay flow, Re=40
   Re  = 40;
   [U,V]=kovasznay(X,Y,Re);
   T=0*U;

%  Taylor vortices
%  U = -sin(pi*X).*cos(pi*Y);
%  V =  cos(pi*X).*sin(pi*Y);

if istep==0; ;       %% Set Desired IC

   U = 0    + 0*X; % U=Mu.*U;
   V = 0    + 0*X; % V=Mv.*V;
   T = 0    + 0*X; % T=Mt.*T;

else;                %% Set Desired BC

   U = Mu.*Uo + (1-Mu).*U;   %% Old value in interior, new value on inflow
   V = Mv.*Vo + (1-Mv).*V;
   T = Mt.*To + (1-Mt).*T;

end;

