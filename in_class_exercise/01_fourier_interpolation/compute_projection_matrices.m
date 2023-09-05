function [F, S, x, y, k] = compute_projection_matrices(N, M, l, r);
% Compute the two operation matrices for Fourier analysis and 
% Fourier synthesis, respectively. Domain x-range is statically
% set a [-pi, pi).
%
% N: Number of sampling points
% M: Number of synthesis points
% l: Left boundary. Defaulted to -pi
% r: Right boundary. Defaulted to pi
%
% Return:
% F: Fourier analysis matrix (N x N)
% M: Fourier synthesis matrix (M x N)
% x: Uniformly-spaced vector spanning l-r (length N)
% y: Uniformly-spaced vector spanning l-r (length M)
% k: Fourier modes (length N)

    if nargin() == 2;
        l = -pi;
        r = pi;
    elseif nargin() == 3;
        r = pi;
    elseif nargin() ~= 4;
        error("Wrong number of arguments");
    end;

    scale_factor = 2 * pi / (r - l);

    x = l + (r - l) / N * [0:N-1];
    y = l + (r - l) / M * [0:M-1];

    k = [-N/2:(N/2-1)] * scale_factor;
    if mod(N, 2) ~= 0;
        k = [-(N-1)/2:(N-1)/2];
    end;

    i = sqrt(-1);
    F = (1 / N) * exp(-i * k' * x);
    S = exp(i * y' * k);
    
end;
