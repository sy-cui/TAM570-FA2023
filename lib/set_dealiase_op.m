function [JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);
    % Set operation matrices for advection dealiasing
    % Output (cell arrays):
    %   JRx :   Interpolated Rx
    %   JD  :   (J * D) Differentiate then interpolate
    %   Bm  :   Mass matrix for finer mesh with Jacobian integrated
    validateattributes(dim,{'numeric'},{'scalar','integer'});
    if  length(D) ~= dim || ...
        length(wm)~= dim || ...
        length(Jm)~= dim || ...
        size(Rx,1)~= dim || ...
        size(Rx,2)~= dim;
        error("Dimension mismatch");
    end;

    if dim == 2; error("Unsupported for now"); end;

    if dim == 3;
        % JRx = diag(Jm) * Rx
        JRx = cell(dim,dim);
        for i=1:dim; for j=1:dim;
            JRx{i,j} = tensor3(Jm{3},Jm{2},Jm{1},Rx{i,j});
        end; end;

        % JD = diag(Jm) * [Dr,Ds,Dt]^T
        JD = cell(1,dim);
        for i=1:dim;
            JD{i} = Jm{i} * D{i};
        end;

        % Bm = (w_i w_j w_k) J * Jac_ijk
        nmx = length(wm{1}); Bmx = reshape(wm{1},nmx,1,1);
        nmy = length(wm{2}); Bmy = reshape(wm{2},1,nmy,1);
        nmz = length(wm{3}); Bmz = reshape(wm{3},1,1,nmz);
        Bm = (Bmx.*Bmy.*Bmz).*tensor3(Jm{3},Jm{2},Jm{1},Jac);
        fprintf("Volume from set_dealiase_op: %.15f\n", sum(sum(sum(Bm))));

    end; % if dim == 3

end;