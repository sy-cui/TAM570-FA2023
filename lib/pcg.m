function [x,k,res,xr] = pcg(rhs,Pinv,S,Linv,R,D,G,B,b0,ndt,max_iter=10000,tol=1e-5);
    k = 0;
    r = rhs;    % Input rhs must be restricted
    x = 0*r;
    z = fdm_solve_3d(r,S,Linv);  % Need to solve
    % z = r;
    p = z;
    rz0 = 1; 
    rz1 = sum(sum(sum(z.*r)));
    res = abs(sqrt(rz1));

    while k < max_iter && res > tol;
        rz0 = rz1;
        Ap = viscous_op(p,R,D,G,B,b0,ndt);
        alpha = rz0 / sum(sum(sum(p.*Ap)));
        x = x + alpha * p;
        r = r - alpha * Ap;
        z = fdm_solve_3d(r,S,Linv);  % Need to solve
        % z = r;
        rz1 = sum(sum(sum(z.*r)));
        res = abs(sqrt(rz1));
        beta = rz1 / rz0;
        p = z + beta * p;
        k = k + 1;
    end;

    if k == 0;
        xr = z;
    else;
        xr = x;
    end;

    x = t3w(R,xr,1);    % Transpose
end;