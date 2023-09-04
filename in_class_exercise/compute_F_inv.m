%% 
function F_inter = compute_F_inv(N, M)
    k = -(N/2):1:(N/2-1);
    temp = linspace(0, N-1, N);
    h = (2*pi)/N;
    xj = -pi + h*(temp);

    kx = k' * xj;
    i = sqrt(-1);

    F = (1/N) * exp(-i * kx);

    %uj = functions(xj);
    %uk = F * (uj');
    hl = (2*pi)/M;
    temp = linspace(0, M-1, M);
    yl = -pi + hl*(temp);

    kyl = yl' * k;

    F_inv = exp(i * kyl);
    F_inter = real(F_inv * F);
end
