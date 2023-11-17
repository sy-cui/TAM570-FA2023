function [X,Y] = gen_box_2d(
    N,          % Legendre order in each element
    Nelx,       % Number of elements in x
    Nely,       % Number of elements in y
    xb=[-1 1],  % x-direciton boundary points
    yb=[-1 1]   % y-direciton boundary points
);
    Le_x = (xb(2) - xb(1)) / Nelx;
    Le_y = (yb(2) - yb(1)) / Nely;
    [z, ~] = zwgll(N);
    zx = 0.5 * (z + 1) * Le_x;
    zy = 0.5 * (z + 1) * Le_y;
    [Zx, Zy] = ndgrid(zx, zy);
    X = zeros(N+1,Nelx*Nely,N+1);
    Y = zeros(N+1,Nelx*Nely,N+1);

    for i=1:Nelx; for j=1:Nely;
        X(:,(i-1)*Nely+j,:) = xb(1) + Zx(:,:) + (i-1)*Le_x;
        Y(:,(i-1)*Nely+j,:) = yb(1) + Zy(:,:) + (j-1)*Le_y;
    end; end;

end;