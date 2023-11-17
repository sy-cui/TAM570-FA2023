function F = tensor3(Az,Ay,Ax,F);

%
%  Note, this code also works if Ax, Ay, or Az are constants instead of matrices
%

[nx,ny,nz]=size(F);

[mx,nxa]=size(Ax);
[my,nya]=size(Ay);
[mz,nza]=size(Az);

if mx==1 && nxa==1; mx=nx; end;    %% Handle the 1x1 Ax case
if my==1 && nya==1; my=ny; end;    %% Handle the 1x1 Ay case
if mz==1 && nza==1; mz=nz; end;    %% Handle the 1x1 Az case

F=reshape(F,nx,ny*nz); F=Ax*F;
F=reshape(F,mx*ny,nz); F=F*Az';

F=reshape(F,mx,ny,mz); 
  if ny==my; for k=1:nz; F(:,:,k) = F(:,:,k)*Ay'; end;
  else; 
    Fy=zeros(mx,my,mz);
    for k=1:mz; Fy(:,:,k) = F(:,:,k)*Ay'; end; F=Fy;
  end;

