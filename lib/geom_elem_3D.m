function [Rx,G,B,Jac]=geom_elem_3D(
    X,Y,Z,...       % Coordinates in x, y, z
    D,...           % Differentiation matrices wrt. r,s,t
    w               % quadrature weights in r,s,t
    );      

    % Compute integration operators for 3D deformed geometry
    % Rx = [dr/dx ds/dx dt/dx
    %       dr/dy ds/dy dt/dy
    %       dr/dz ds/dz dt/dz]
    % Xr = [dx/dr dx/ds dx/dt
    %       dy/dr dy/ds dy/dt
    %       dz/dr dz/ds dz/dt]
    % Xr * Rx' = I
    % Jac = det(Xr) = 1/det(Rx) = 
    %   dx/dr * (dy/ds * dz/dt - dz/ds * dy/dt) -
    %   dx/ds * (dy/dr * dz/dt - dz/dr * dy/dt) +
    %   dx/dt * (dy/dr * dz/ds - dz/dr * dy/ds)
    % G = Rx' * JB * Rx
    %
    % Output:
    %   Rx  :   3x3 cell array containing dr_j / dx_i
    %   G   :   G = Rx^T diag(B) Rx. 1x6 cell array containing
    %           { Grr, Gss, Gtt, Grs, Grt, Gst }
    %   B   :   Mass matrix integrated with Jacobian
    %   Jac :   Volumetric Jacobian 

    [wr,ws,wt] = w{:};
    validateattributes(wr, {'numeric'}, {'vector'});
    validateattributes(ws, {'numeric'}, {'vector'});
    validateattributes(wt, {'numeric'}, {'vector'});

    nx = length(wr); ny = length(ws); nz = length(wt);
    validateattributes(X, {'numeric'}, {'3d', 'size', [nx, ny, nz]});
    validateattributes(Y, {'numeric'}, {'3d', 'size', [nx, ny, nz]});
    validateattributes(Z, {'numeric'}, {'3d', 'size', [nx, ny, nz]});

    [Dhr,Dhs,Dht] = D{:};
    validateattributes(Dhr, {'numeric'}, {'2d', 'size', [nx, nx]});
    validateattributes(Dhs, {'numeric'}, {'2d', 'size', [ny, ny]});
    validateattributes(Dht, {'numeric'}, {'2d', 'size', [nz, nz]});

    % Compute dx_i / dr_j
    dxdr=tensor3(1,1,Dhr,X); dxds=tensor3(1,Dhs,1,X); dxdt=tensor3(Dht,1,1,X);
    dydr=tensor3(1,1,Dhr,Y); dyds=tensor3(1,Dhs,1,Y); dydt=tensor3(Dht,1,1,Y);
    dzdr=tensor3(1,1,Dhr,Z); dzds=tensor3(1,Dhs,1,Z); dzdt=tensor3(Dht,1,1,Z);

    % Compute Jacobian
    Jac = (
          dxdr.*(dyds.*dzdt - dzds.*dydt)
        - dxds.*(dydr.*dzdt - dzdr.*dydt)
        + dxdt.*(dydr.*dzds - dzdr.*dyds)
    );
    Ji  = 1./Jac;
    Ji(abs(Ji) > 1e14) = 0; % Remove singularities

    % Compute dr_j / dx_i
    drdx= Ji.*(dyds.*dzdt-dzds.*dydt); dsdx=-Ji.*(dydr.*dzdt-dzdr.*dydt); dtdx= Ji.*(dydr.*dzds-dzdr.*dyds); 
    drdy=-Ji.*(dxds.*dzdt-dzds.*dxdt); dsdy= Ji.*(dxdr.*dzdt-dzdr.*dxdt); dtdy=-Ji.*(dxdr.*dzds-dzdr.*dxds); 
    drdz= Ji.*(dxds.*dydt-dyds.*dxdt); dsdz=-Ji.*(dxdr.*dydt-dydr.*dxdt); dtdz= Ji.*(dxdr.*dyds-dydr.*dxds); 

    % Compute diagonal mass matrix
    Bx=reshape(wr,nx,1,1); By=reshape(ws,1,ny,1); Bz=reshape(wt,1,1,nz);
    B = Jac.*(Bx.*By.*Bz);
    fprintf("Volume from geom_elem_3D: %.15f\n", sum(sum(sum(B))));

    Grr=B.*(drdx.*drdx+drdy.*drdy+drdz.*drdz);
    Grs=B.*(dsdx.*drdx+dsdy.*drdy+dsdz.*drdz); Gss=B.*(dsdx.*dsdx+dsdy.*dsdy+dsdz.*dsdz);
    Grt=B.*(dtdx.*drdx+dtdy.*drdy+dtdz.*drdz); Gst=B.*(dtdx.*dsdx+dtdy.*dsdy+dtdz.*dsdz); Gtt=B.*(dtdx.*dtdx+dtdy.*dtdy+dtdz.*dtdz);

    Rx = {  drdx dsdx dtdx;
            drdy dsdy dtdy;
            drdz dsdz dtdz  };
    G = {Grr,Gss,Gtt,Grs,Grt,Gst};
end;

