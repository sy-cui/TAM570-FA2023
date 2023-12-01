function [dA] = diag_a_2d(Dr,Ds,Grr,Gss,Grs);
    % Return the diagonal of the Neumann operator A in grid form
    dA = (Dr.*Dr)'*Grr + Gss*(Ds.*Ds);
    dA = dA + 2 * diag(diag(Dr)) * Grs * diag(diag(Ds));
end;


