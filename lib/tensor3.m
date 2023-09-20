function F = tensor3(Az,Ay,Ax,F);

%
%  Note, this code also works if Ax, Ay, or Az are constants instead of matrices
%

[nx,ny,nz]=size(F);

F=reshape(F,nx,ny*nz); F=Ax*F;
F=reshape(F,nx*ny,nz); F=F*Az';
F=reshape(F,nx,ny,nz); for k=1:nz; F(:,:,k) = F(:,:,k)*Ay'; end;



