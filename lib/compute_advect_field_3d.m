function [Cr,Cs,Ct] = compute_advect_field_3d(
    U,V,W,...   % Velocities in x, y, z
    Jm,...      % Interpolation matrices to fine mesh
    JRx         % Interpolated Rx
    );
    % Compute advecting velociy field in r,s,t
    JU = tensor3(Jm{3},Jm{2},Jm{1},U);
    JV = tensor3(Jm{3},Jm{2},Jm{1},V);
    JW = tensor3(Jm{3},Jm{2},Jm{1},W);
    Cr = JU.*JRx{1,1} + JV.*JRx{2,1} + JW.*JRx{3,1}; % Cx*dr/dx + Cy*dr/dy + Cz*dr/dz 
    Cs = JU.*JRx{1,2} + JV.*JRx{2,2} + JW.*JRx{3,2}; % Cx*dr/dx + Cy*dr/dy + Cz*dr/dz 
    Ct = JU.*JRx{1,3} + JV.*JRx{2,3} + JW.*JRx{3,3}; % Cx*dr/dx + Cy*dr/dy + Cz*dr/dz 

end;