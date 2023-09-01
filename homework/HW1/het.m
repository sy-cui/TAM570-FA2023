function [T, dT] = het(r, x, k, p, L, bessel_roots);
% Compute the mode at N
% r: radial coordinates
% x: axial coordiantes
% p: quadrature order 
% L: cylinder length
% k: mode index
% bessel_roots: roots of J_0

    [z, w] = zwgll(p); 
    z = 0.5 * (z + 1); w = 0.5 * w; % Scale to [0, 1]

    j0 = besselj(0, bessel_roots(k) .*z);
    sh_by_ch = tanh(bessel_roots(k) * L);
    beta_k = sum(w.*j0.*z) / sum(w.*j0.^2.*z);

    T = (
        beta_k 
        / bessel_roots(k) 
        * (sh_by_ch * cosh(bessel_roots(k) * x) - sinh(bessel_roots(k) * x)) 
        * besselj(0, bessel_roots(k)*r)
    );
    dT = (
        beta_k 
        * (sh_by_ch * sinh(bessel_roots(k) * x) - cosh(bessel_roots(k) * x)) 
        * besselj(0, bessel_roots(k)*r)
    );

end;