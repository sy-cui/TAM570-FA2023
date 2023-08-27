function [xb, ub, us] = ad1d(n, Pe, c, L);
    % n :   number of grid points
    % Pe:   Peclet number (advection / diffusion)
    % c :   Advective speed
    % L :   Domain length
    
    nu = c * L / Pe;    % Diffusion coefficient
    h = L/n;            % Grid spacing

    n1=n-1; e = ones(n1,1);
    A = -nu*spdiags([e -2*e e], -1:1, n1, n1)/(h*h);    % A - SPD
    C = c*spdiags([-e 0*e e], -1:1, n1, n1)/(2*h);      % C - skew-symm
    H = A+C; % Advection-Diffusion operator
    x = [1:n1]*h; x=x';
    f = 1+0*x;
    u = H\f;
    ub = [0; u; 0]; xb=[0; x; L]; % Extend u and x to boundaries
    xi = xb / L;
    us = L / c * (xi - (exp(Pe*(xi - 1)) - exp(-Pe))/(1 - exp(-Pe)));   % Analytical solution