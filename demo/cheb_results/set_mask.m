function [mask]=set_mask(bc_left,bc_right,bc_lower,bc_upper,N,Nelx,Nely)

N1=N+1;
E =Nelx*Nely;

% disp([N E Nelx Nely])

mask = ones(N1,Nelx,Nely,N1);

if bc_left =='D'; mask( 1,   1,:,:)=0; end;
if bc_right=='D'; mask(N1,Nelx,:,:)=0; end;
if bc_lower=='D'; mask(:,:,   1, 1)=0; end;
if bc_upper=='D'; mask(:,:,Nely,N1)=0; end;


mask=reshape(mask,N1,E,N1);
