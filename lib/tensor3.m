function G = tensor3(Az,Ay,Ax,F);

%
%  Note, this code also works if Ax, Ay, or Az are constants instead of matrices
%

[nx,ny,nz]=size(F);

F=reshape(F,nx,ny*nz); F=Ax*F; nx=size(F)(1);
F=reshape(F,nx*ny,nz); F=F*Az'; nz=size(F)(end);
G = zeros(nx,size(Ay)(1),nz);

F=reshape(F,nx,ny,nz); 
for k=1:nz; G(:,:,k) = F(:,:,k)*Ay'; end;


