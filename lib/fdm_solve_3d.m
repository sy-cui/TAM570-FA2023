function [G] = fdm_solve_3d(F,S,Linv);
    % Fast Diagonalization Method Poisson solve

    % G = tensor3(S{3},S{2},S{1},Linv.*tensor3(S{3}',S{2}',S{1}',F));
    G = t3w(S,Linv.*t3w(S,F,1));

end;