function [Gx,Gy,Gz] = grad_3d(F,Rx,D);
    % Compute gradient of scalar field F
    Fr = tensor3(1,1,D{1},F); % df/dr
    Fs = tensor3(1,D{2},1,F); % df/ds
    Ft = tensor3(D{3},1,1,F); % df/ds

    Gx = Fr.*Rx{1,1} + Fs.*Rx{1,2} + Ft.*Rx{1,3};
    Gy = Fr.*Rx{2,1} + Fs.*Rx{2,2} + Ft.*Rx{2,3};
    Gz = Fr.*Rx{3,1} + Fs.*Rx{3,2} + Ft.*Rx{3,3};

end;