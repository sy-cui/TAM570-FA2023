function [Xout, Yout] = laplace_blend_2d(Xin, Yin);
    % Xin, Yin must contain the boundary points
    if size(Xin) ~= size(Yin);
        error("Xin, Yin dimension mismatch");
    end;

    [nx, ny] = size(Xin);
    [Ahx, Bhx, Chx, Dhx, zx, wx] = semhat(nx-1);
    [Ahy, Bhy, Chy, Dhy, zy, wy] = semhat(ny-1);
    Ihx = speye(nx); Ihy = speye(ny);
    Rx = Ihx(2:end-1, :); Ry = Ihy(2:end-1, :);

    Ax = Rx * Ahx * Rx'; Ay = Ry * Ahy * Ry';
    Bx = Rx * Bhx * Rx'; By = Ry * Bhy * Ry';

    mask = Xin * 0;
    mask(1,:)=1; mask(:,1)=1; mask(end,:)=1; mask(:,end)=1;
    Mx = mask .* Xin; My = mask .* Yin;
    RHSx = -(Ahx*Mx*Bhy' + Bhx*Mx*Ahy');
    RHSy = -(Ahx*My*Bhy' + Bhx*My*Ahy');

    [Sx,Lx]=gen_eig_decomp(Ax,Bx); Sx = Rx'*Sx;
    [Sy,Ly]=gen_eig_decomp(Ay,By); Sy = Ry'*Sy;
    Linv = 1./(full(diag(Lx))+full(diag(Ly))');

    Xout = Mx + tensor2(Sy, Sx, Linv.*tensor2(Sy', Sx', RHSx));
    Yout = My + tensor2(Sy, Sx, Linv.*tensor2(Sy', Sx', RHSy));
 

