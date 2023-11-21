function [G] = assemble(F,N,Nelx,Nely);
    [nx, ne, ny] = size(F);
    if nx~=N+1 || ny~=N+1 || ne~=Nelx*Nely;
        error("Incompatible shape!")
    end;

    nx = Nelx*N+1;
    ny = Nely*N+1;
    G = zeros(nx,ny);
    for j=1:Nely; for i=1:Nelx;
        el = (j-1)*Nelx+i;
        G((i-1)*N+1:i*N+1, (j-1)*N+1:j*N+1) = squeeze(F(:,el,:));
    end; end;

end;