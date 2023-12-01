function [x,Xk,k,res] = pcg_prj(
    rhs,            % Unrestricted RHS to be solved (b^n)
    Xk,             % Collection of A-conjugate basis vectors. May be empty
    Pinv,           % Inverse of preconditioner matrix P
    S,
    Linv,
    R,              % Restriction matrices. May contain periodic
    D,              % Differentiation matrices
    G,              % G matrix = Rx^T B Rx
    B,              % Unrestricted mass matrix
    b0,             % First BDFk parameter
    ndt,            % nu * dt
    iter,           % Time loop iterations.
    sdim=10,        % Subspace dimension
    max_iter=10000, % Maximum PCG iteration
    tol=1e-5        % PCG tolerance
    );
    % A-conjugate projection method for solving successive linear systems 
    % using preconditioned conjugate gradient
    % See https://www.sciencedirect.com/science/article/pii/S0045782598000127

    if isempty(Xk) || iter==1; % Solve directly
        [x,k,res,xr]=pcg(t3w(R,rhs),Pinv,S,Linv,R,D,G,B,b0,ndt,max_iter,tol);
        Axr = viscous_op(xr,R,D,G,B,b0,ndt);
        xr_anorm = sqrt(sum(sum(sum(xr.*Axr))));
        if xr_anorm > 1e-14; Xk = [reshape(xr,[],1)/xr_anorm]; end;
        return;
    end;

    R_rhs = t3w(R,rhs);             % Restricted RHS
    b = reshape(R_rhs,[],1);        % Restricted RHS as a column vector
    alpha = Xk'*b;                  % (xt, b) coefficients

    % \bar{x} = sum_j { (xt_j, b)*xt_j }
    xb = reshape(Xk*(Xk'*b),size(R_rhs)); 

    % f = b - A\bar{x} (restricted)                    
    f = R_rhs-viscous_op(xb,R,D,G,B,b0,ndt);

    % Solve to tol. Get restricted soln only
    [~,k,res,y]=pcg(f,Pinv,S,Linv,R,D,G,B,b0,ndt,max_iter,tol);

    % x = y+xb prolongated          
    xr = y+xb; x = t3w(R,xr,1);

    if mod(iter, sdim) == 1; % reset basis
        Ax = viscous_op(xr,R,D,G,B,b0,ndt);
        x_anorm = sqrt(sum(sum(sum(xr.*Ax))));
        if x_anorm > 1e-14; Xk = [reshape(xr,[],1)/x_anorm]; end;

    else;
        % \hat{f} = Ay as column vector                   
        fh = reshape(viscous_op(y,R,D,G,B,b0,ndt),[],1);

        % Reshape y to colume vector
        y = reshape(y, [], 1);

        % beta_j = (xt_j, fh)
        beta = Xk' * fh;

        % \hat{x} = y - sum_j {beta_j * xt_j}
        xh = y - Xk*beta;        

        % A-norm of xh: sqrt{ (y, Ay) - sum_j {beta_j^2} }
        xh_anorm = sqrt(y'*fh - sum(beta.^2));  

        % Append xh / ||xh||_A to basis set
        if xh_anorm > 1e-14; Xk = [Xk xh/xh_anorm]; end;
        
    end; % endif

end; % pcg_prj