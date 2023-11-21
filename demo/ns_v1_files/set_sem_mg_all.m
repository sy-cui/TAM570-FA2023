function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,...
          Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA,...
          Dc,Jc,Grrc,Grsc,Gssc,Blc,Qc,glo_numc,Mc,dAc] = set_sem_mg_all(N);



Nc=ceil(N/2)+1;
Nc=N-1;

[U,V,T,zc,wc,Dc,X,Y,Grrc,Grsc,Gssc,Blc,Xr,Rx,Jac,Qc,glo_numc,Mu,Mv,Mc,Mt,ifnull,unxa_v,unya_v,BC_all,dAc]...
             = set_sem_all(Nc);

[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]...
             = set_sem_all(N);

Jc=interp_mat(z,zc);


