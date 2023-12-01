function [R] = set_restriction(N,bcs);
    % Return the 1D restriction matrices for corresponding boundary condition
    % bcs must be a cell array of characters.
    % 'D' or 'd' represent Dirichlet
    % 'P' or 'p' represent periodic
    % Anything else will be parse as Neumann
    validateattributes(N,{'numeric'},{'vector','integer','positive'});
    dim = length(N);
    if length(bcs) ~= 2*dim;
        error("Incompatible numbers of boundary conditions");
    end;
    bcs = lower(bcs);

    R = cell(1,dim);
    for i=1:dim;
        if bcs{2*i-1} == 'p';
            if bcs{2*i} ~= 'p'; error('Inconsistent periodic bc'); end; 
            R{i} = r_periodic(N(i));
            continue;
        end;

        R{i} = speye(N(i)+1);
        if bcs{2*i} == 'p'; error('Inconsistent periodic bc'); end; 
        if bcs{2*i-1} == 'd'; R{i} = R{i}(2:end,:); end;
        if bcs{2*i} == 'd'; R{i} = R{i}(1:end-1,:); end;

    end;
end;