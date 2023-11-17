function Wb = axb(Ub,b0,nu,Bb,Grr,Grs,Gss,Dhr,Dhs);

%
%   A = D^T G D
%

% Ub  - rhs, in grid form (full grid)
% Gij - geometric factors, in grid form
% Dhr - derivative in Omega-hat, r direction
% Dhs - derivative in Omega-hat, s direction


Ur = Dhr*Ub;
Us = Ub*Dhs';

Wr = nu.*(Grr.*Ur + Grs.*Us);  % Can support variable nu
Ws = nu.*(Grs.*Ur + Gss.*Us);
Wb = Dhr'*Wr + Ws*Dhs + b0.*(Bb.*Ub);

