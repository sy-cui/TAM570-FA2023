function [U,V,T] = set_ic_bc(Uo,Vo,To,Mu,Mv,Mt,X,Y,IC);

%  Kovasznay flow, Re=40
   Re  = 40;
   [U,V]=kovasznay(X,Y,Re);

%  Taylor vortices
%  U = -sin(pi*X).*cos(pi*Y);
%  V =  cos(pi*X).*sin(pi*Y);

if IC==1; ;          %% Set Desired IC
   U = 0 + 0*X;
   V = 0 + 0*X;
   T = 0 + 0*X;
else;                %% Set Desired BC
%  U = 1 - Y.*Y;  %% Desired field at inflow
%  V = 0 + 0*X;   %% Desired field at inflow
   T = 1 + 0*X;   %% Desired field at inflow

   U = Mu.*Uo + (1-Mu).*U;   %% Old value in interior, new value on inflow
   V = Mv.*Vo + (1-Mv).*V;
   T = Mt.*To + (1-Mt).*T;
end;

