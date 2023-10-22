hdr

N=5; n=N+1; 

M=ceil(21.2*N); [zf,w]=zwuni(M); 

[zu,w]=zwuni(N); 
[zg,w]=zwgll(N);

Ju=f_interp_mat(zf,zu);
Jg=  interp_mat(zf,zg);

hold off;
z=zg; J=Jg;
for j=1:n
    s=1.5*j;
    f=zeros(n,1)+s; f(j)=1+s; ff=J*f;
    plot(2*z,0*z+s,'k.-',ms,14,...
         2*z,f    ,'r.',ms,14,...
         2*z,0*z+s,'k.-',ms,14,...
         2*zf,ff  ,'r-',lw,2); hold on;
end;
axis equal; axis off; 
print -dpng 'legendre_x5.png'

hold off;
z=zu; J=Ju;
for j=1:n
    s=1.5*j;
    f=zeros(n,1)+s; f(j)=1+s; ff=J*f;
    fp=f;
    if j==1 fp(1)=0.5+s; fp(n)=0.5+s; end;
    if j==n fp(1)=0.5+s; fp(n)=0.5+s; end;
    plot(2*z,0*z+s,'k.-',ms,14,...
         2*z,fp   ,'r.',ms,14,...
         2*z,0*z+s,'k.-',ms,14,...
         2*zf,ff  ,'r-',lw,2); hold on;
end;
axis equal; axis off; 
print -dpng 'fourier_x5.png'
