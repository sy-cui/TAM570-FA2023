function [Xs,Ys,Zs] = morph_sphere(X,Y,Z,ri,ro);
    % Morph domain into a sphere of inner radius ri
    % and outer radius ro centered at (0,0)

    validateattributes(ri,{'numeric'},{'scalar','positive'});
    validateattributes(ro,{'numeric'},{'scalar','positive'});

    if ri >= ro; 
        error("Inner radius must be smaller than outer radius");
    end;

    x0=min(min(min(X))); x1=max(max(max(X)));  Lx=x1-x0;
    y0=min(min(min(Y))); y1=max(max(max(Y)));  Ly=y1-y0;
    z0=min(min(min(Z))); z1=max(max(max(Z)));  Lz=z1-z0;

    X = (2.0 / Lx) * (X - x0) - 1.0;    % On [-1,1]
    Y = (1.0 / Ly) * (Y - y0);          % On [0, 1]
    Z = (1.0 / Lz) * (Z - z0);          % On [0, 1]

    azi = pi*X;                         % Azimuthal
    pol = pi*(1-Y);                     % Polar
    rad = ri + (ro - ri) * Z;           % Radial

    Xs=rad.*sin(pol).*cos(azi);
    Ys=rad.*sin(pol).*sin(azi);
    Zs=rad.*cos(pol);

end;
