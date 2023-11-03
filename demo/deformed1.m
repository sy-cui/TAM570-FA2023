addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

hold off

N=9;
Nf=100;
N2=2;

[z,w]=zwgll(N);
[zf,wf]=zwgll(Nf);
[z2,w2]=zwgll(N2);

Jf=interp_mat(zf,z2);
J =interp_mat(z ,z2);

[X2,Y2]=ndgrid(z2,z2);

Y2(1,1) = -1.1;
Y2(2,1) = -0.95;
Y2(2,3) =  1.05;

X2(1,1) = -0.9;
X2(2,1) =  0.05;

X2(1,2) = -1.25;
X2(3,2) =  0.80;

Y2=Y2+1;

Xbar=sum(sum(X2))/9; X2(2,2) = Xbar;
Ybar=sum(sum(Y2))/9; Y2(2,2) = Ybar;

X = J*X2*J'; Y = J*Y2*J';
Xf= Jf*X2*J'; Yf = Jf*Y2*J';
Xg= J*X2*Jf'; Yg = J*Y2*Jf';

plot(X,Y,'ko',Xf,Yf,'r-',Xg',Yg','b-',lw,1.5);

text(X2(1,2)-.20,Y2(1,2)-.00,'$a$',intp,ltx,fs,28)
text(X2(3,2)+.25,Y2(3,2)+.00,'$b$',intp,ltx,fs,28)
text(X2(2,1)-.05,Y2(2,1)-.25,'$c$',intp,ltx,fs,28)
text(X2(2,3)+.00,Y2(2,3)+.20,'$d$',intp,ltx,fs,28)

dX=max(max(X))-min(min(X)); dX=1.25*dX/2;
dY=max(max(Y))-min(min(Y)); dY=1.35*dY/2; dX=max(dX,dY);
aX=sum(sum(X))/sum(sum(0*X+1));
aY=sum(sum(Y))/sum(sum(0*Y+1));
axis([aX-dX-.1 aX+dX-.1 aY-dX aY+dX])

axis square;
pause
% print -dpng deformed1.png

