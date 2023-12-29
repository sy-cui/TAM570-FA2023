function se_quiver(X,Y,U,V,name);
hdr

E=size(X,2);
N1=size(X,1);

% hold off;
for e=1:E
    x=X(:,e,:); x=reshape(x,N1,N1);
    y=Y(:,e,:); y=reshape(y,N1,N1);
    u=U(:,e,:); u=reshape(u,N1,N1);
    v=V(:,e,:); v=reshape(v,N1,N1);
    u2=sqrt(max(max(u.*u+v.*v)));
    h=quiver(x,y,u,v,'autoscalefactor',u2); 
    hold on;
end;
title(name,fs,16);
xlabel('X',fs,16);
ylabel('Y',fs,16);
% view(55,33);
