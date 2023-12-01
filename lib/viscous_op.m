function [HF] = viscous_op(F,R,D,G,B,b0,ndt);
    % Compute LHS viscous operation
    % Assume F is restricted.
    % To apply to unrestricted F, use R = {1,1,1}.
    % R (b0*B + nu*dt*D^T G D) R^T F
    %   G{1} = Grr
    %   G{2} = Gss
    %   G{3} = Gtt
    %   G{4} = Grs = Gsr
    %   G{5} = Grt = Gtr
    %   G{6} = Gst = Gts

    % Prolongate
    F = t3w(R,F,1);

    %% D * F
    Fr = tensor3(1, 1, D{1}, F); 
    Fs = tensor3(1, D{2}, 1, F); 
    Ft = tensor3(D{3}, 1, 1, F); 

    %% G * D * F
    GFr = G{1}.*Fr + G{4}.*Fs + G{5}.*Ft;
    GFs = G{4}.*Fr + G{2}.*Fs + G{6}.*Ft;
    GFt = G{5}.*Fr + G{6}.*Fs + G{3}.*Ft;

    %% (b0 * B + ndt * D^t * G * D) * F
    HF = b0*B.*F + ndt.*(...
          tensor3(1,1,D{1}',GFr)...
        + tensor3(1,D{2}',1,GFs)...
        + tensor3(D{3}',1,1,GFt));

    %% Restrict
    HF = t3w(R,HF);

end;