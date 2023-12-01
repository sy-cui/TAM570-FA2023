function [rhs] = pp_rhs(Ut,Vt,Wt,Rx,D,dt);
    % Compute RHS of the pressure Poisson equation
    % \int_\Omega Grad{q} * Ut dV * (1/dt)
    % This is equivalent to
    % q^T D^T Rx^T diag(B) ut / dt
    % B is assumed to be pre-applied to Ut

    Ur = Ut.*Rx{1,1} + Vt.*Rx{2,1} + Wt.*Rx{3,1};
    Us = Ut.*Rx{1,2} + Vt.*Rx{2,2} + Wt.*Rx{3,2};
    Ut = Ut.*Rx{1,3} + Vt.*Rx{2,3} + Wt.*Rx{3,3};

    rhs = (1/dt)*(
          tensor3(1,1,D{1}',Ur)...
        + tensor3(1,D{2}',1,Us)...
        + tensor3(D{3}',1,1,Ut)
    );

end;