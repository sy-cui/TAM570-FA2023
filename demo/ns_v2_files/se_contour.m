function se_contour(X,Y,Z,name,ncontour);
hdr

E=size(X,2);
N1=size(X,1);

zm=min(min(min(Z)));
zM=max(max(max(Z)));

% ncontour = 20;
zc=zm + (zM-zm)*[0:ncontour]'/ncontour;

% hold off;
for e=1:E
    x=X(:,e,:); x=reshape(x,N1,N1);
    y=Y(:,e,:); y=reshape(y,N1,N1);
    z=Z(:,e,:); z=reshape(z,N1,N1);
    contour(x,y,z,zc,lw,2); hold on;
end;

name = [name ' Min/Max: ' num2str([zm zM]) ];
title(name,fs,16);
xlabel('X',fs,16);
ylabel('Y',fs,16);
