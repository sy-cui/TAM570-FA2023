function [Grr,Grs,Gss,B,Jac]=geom_elem(X,Y,Dh,w);
    % Compute integration operators for deformed geometry
    % Assume for each element, Nx = Ny = N
    % Rx = [dr/dx ds/dx
    %       dr/dy ds/dy]
    % Xr = [dx/dr dx/ds
    %       dy/dr dy/ds]
    % Xr * Rx' = I
    % Jac = det(Xr) = 1/det(Rx) = dx/dr * dy/ds - dy/dr * dx/ds
    % G = Rx' * JB * Rx

    n = length(w);  N=n-1;
    NE = size(X,2);

    validateattributes(X, {'numeric'}, {'3d', 'size', [n, NE, n]});
    validateattributes(Y, {'numeric'}, {'3d', 'size', [n, NE, n]});
    validateattributes(Dh, {'numeric'}, {'2d', 'size', [n, n]});

    Be=w*w';
    for e=1:NE
        B(:,e,:)=Be(:,:);
    end;

    %  Compute dx_i / dr_j
    dxdr = tensor3(1,1,Dh,X);
    dxds = tensor3(Dh,1,1,X);
    dydr = tensor3(1,1,Dh,Y);
    dyds = tensor3(Dh,1,1,Y);

    Jac = dxdr.*dyds - dxds.*dydr;
    Ji  = 1./Jac;

    drdx =  Ji.*dyds;
    drdy = -Ji.*dxds;
    dsdx = -Ji.*dydr;
    dsdy =  Ji.*dxdr;

    B = Jac.*B; % nr x ns diagonal mass matrix on Omega-hat, incl. Jac

    Grr=B.*(drdx.*drdx + drdy.*drdy);
    Grs=B.*(drdx.*dsdx + drdy.*dsdy);
    Gss=B.*(dsdx.*dsdx + dsdy.*dsdy);

end;

