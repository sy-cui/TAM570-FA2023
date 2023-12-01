function [Xc,Yc] = morph_circle(X,Y,ri,ro);
    % Morph domain into a circle of inner radius ri
    % and outer radius ro centered at (0,0)

    validateattributes(ri,{'numeric'},{'scalar','positive'});
    validateattributes(ro,{'numeric'},{'scalar','positive'});

    if ri >= ro; 
        error("Inner radius must be smaller than outer radius");
    end;

    x0=min(min(min(X))); x1=max(max(max(X)));  Lx=x1-x0;
    y0=min(min(min(Y))); y1=max(max(max(Y)));  Ly=y1-y0;

    X = (2.0 / Lx) * (X - x0) - 1.0;    % On [-1,1]
    Y = (1.0 / Ly) * (Y - y0);          % On [0, 1]

    dr = ro - ri;

    Theta = pi*X;
    R     = ri + dr*Y;

    Xc=R.*sin(Theta);
    Yc=R.*cos(Theta);

end;
