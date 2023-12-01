function [Cx,Cy,Cz]=curl_3d(U,V,W,Rx,D);
    % Compute the curl of vector field [U,V,W]
    %   Cx = dw/dy - dv/dz
    %   Cy = du/dz - dw/dx
    %   Cz = dv/dx - du/dy
    %   Rx = {  drdx dsdx dtdx 
    %           drdy dsdy dtdy
    %           drdz dsdz dtdz  }
    % Output: Cx, Cy, Cz

    dudr=tensor3(1,1,D{1},U); duds=tensor3(1,D{2},1,U); dudt=tensor3(D{3},1,1,U);
    dvdr=tensor3(1,1,D{1},V); dvds=tensor3(1,D{2},1,V); dvdt=tensor3(D{3},1,1,V);
    dwdr=tensor3(1,1,D{1},W); dwds=tensor3(1,D{2},1,W); dwdt=tensor3(D{3},1,1,W);

    Cx = (dwdr.*Rx{2,1} + dwds.*Rx{2,2} + dwdt.*Rx{2,3}) - ...  % dw/dy
         (dvdr.*Rx{3,1} + dvds.*Rx{3,2} + dvdt.*Rx{3,3});       % dv/dz
    Cy = (dudr.*Rx{3,1} + duds.*Rx{3,2} + dudt.*Rx{3,3}) - ...  % du/dz
         (dwdr.*Rx{1,1} + dwds.*Rx{1,2} + dwdt.*Rx{1,3});       % dw/dx
    Cz = (dvdr.*Rx{1,1} + dvds.*Rx{1,2} + dvdt.*Rx{1,3}) - ...  % dv/dx
         (dudr.*Rx{2,1} + duds.*Rx{2,2} + dudt.*Rx{2,3});       % du/dy

end;