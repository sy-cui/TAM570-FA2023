function [J] = new_interp_mat(xo, xi);

    nx = length(xi);
    ny = length(xo);

    x = reshape(xi, nx, 1);
    y = reshape(xo, ny, 1);

    xi_m_xj = (x - x') + eye(nx);
    yi_m_xj = y - x';

    alpha_i = prod(xi_m_xj, 2);
    alpha_i_inv = 1 ./ alpha_i';

    J = zeros(ny, nx);

    for i = 1:nx;
        yi_m_xj_temp = yi_m_xj;
        yi_m_xj_temp(:, i) = 1;
        J(:, i) = prod(yi_m_xj_temp, 2);
    end;

    J = J.*alpha_i_inv;

