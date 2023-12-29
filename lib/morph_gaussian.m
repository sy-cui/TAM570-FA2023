function [X,Y,Z] = morph_gaussian(X,Y,Z,a,s,cx,cy);
    % x0=min(min(min(X))); x1=max(max(max(X)));  Lx=x1-x0;
    % y0=min(min(min(Y))); y1=max(max(max(Y)));  Ly=y1-y0;
    z0=min(min(min(Z))); z1=max(max(max(Z)));  Lz=z1-z0;

    % Map to [-1 1]
    % X = (2 / Lx) * (X - x0) - 1;
    % Y = (2 / Ly) * (Y - y0) - 1;
    Z = (2 / Lz) * (Z - z0) - 1;

    R2 = (X-cx).^2 + (Y-cy).^2; 
    s2_inv = 1 / s^2;
    Z = Z + a*0.5*(1-Z).*exp(-R2*s2_inv);

    % Xg = (Lx / 2) * (X + 1) + x0;
    % Yg = (Ly / 2) * (Y + 1) + y0;
    Z = (Lz / 2) * (Z + 1) + z0;

end;