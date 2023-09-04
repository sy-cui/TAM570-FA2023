%% 
function F_inter = compute_F_inv(N, M)
    if mod(N, 2) == 0
        temp1 = linspace(-N/2, N/2-1, N);
        temp2 = linspace(-M/2, M/2-1, M); 
    else
        temp1 = linspace(-(N-1)/2, (N-1)/2, N);
        temp2 = linspace(-M/2, M/2-1, M); 
    end
    k = -(N/2):1:(N/2-1);
    h = (2*pi)/N;
    xj = h*(temp1);

    kx = k' * xj;
    i = sqrt(-1);

    F = (1/N) * exp(-i * kx);

    %uj = functions(xj);
    %uk = F * (uj');
    hl = (2*pi)/M;
    yl = hl*(temp2);

    kyl = yl' * k;

    F_inv = exp(i * kyl);
    F_inter = real(F_inv * F);
end
