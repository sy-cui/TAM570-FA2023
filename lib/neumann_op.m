function [A,S,L] = neumann_op(
    dim,% dimension of the problem, 2 or 3
    R,  % {Rr, Rs, Rt} restriction matrices
    w,  % {wr, ws, wt} quadrature weights
    D,  % {Dr, Ds, Dt} differentiation matrices
    dl  % {Lr, Ls, Lt} domain lengths
    );
    % Return the restricted Neumann opeartor without deformation
    %       A = R D^T B D R^T
    % in \hat{\Omega} space
    % Also return the generalized eigen-decomposition

    if dim ~= 2 && dim ~= 3;
        error("dim: Dimension must be 2 or 3");
    end;
    
    A = cell(1,dim);
    S = cell(1,dim);
    L = cell(1,dim);

    for i=1:dim;
        Bh = (dl{i}/2)*diag(w{i}); Ah = (4/dl{i}^2)*D{i}'*Bh*D{i};
        A{i} = R{i}*Ah*R{i}';
        [S{i},L{i}] = gen_eig_decomp(A{i},R{i}*Bh*R{i}');
        L{i} = diag(full(L{i}));
    end;
    if dim == 2
        L = reshape(L{1},[],1)+reshape(L{2},1,[]);
    else;
        L = reshape(L{1},[],1,1)+reshape(L{2},1,[],1)+reshape(L{3},1,1,[]);
    end


end;