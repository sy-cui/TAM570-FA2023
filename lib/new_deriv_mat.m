function [D] = new_deriv_mat(x);

    n = length(x);
    xi = reshape(x, n, 1);
    diff_ij = (xi - xi') + eye(n);

    diff_ij_inv = 1 ./ diff_ij;
    D_diag_m1 = sum(diff_ij_inv, 2) - 2;

    alpha_i = prod(diff_ij, 2);   
    alpha_i_inv = 1 ./ alpha_i; 

    D = (alpha_i * alpha_i_inv') ./ diff_ij + diag(D_diag_m1);
