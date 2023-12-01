function [adv] = advect_3d(F,Cr,Cs,Ct,Jm,JD,Bm);
    % Advect the field F with dealiased advecting field
    % Cr, Cs, Ct are outputs from function
    % compute_advect_field_3d
    Dur = tensor3(Jm{3},Jm{2},JD{1},F); % du/dr
    Dus = tensor3(Jm{3},JD{2},Jm{1},F); % du/ds
    Dut = tensor3(JD{3},Jm{2},Jm{1},F); % du/dt
    cdU = Cr.*Dur + Cs.*Dus + Ct.*Dut;
    adv = tensor3(Jm{3}',Jm{2}',Jm{1}',Bm.*cdU);
end;