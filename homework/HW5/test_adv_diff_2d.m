addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

x0 = 0.5; y0 = 0;
nx = 80;
ny = 80;
nu = 0.01;
fcx = @(x,y) -cos(0.5*pi*x).*sin(0.5*pi*y);
fcy = @(x,y) sin(0.5*pi*x).*cos(0.5*pi*y);
f0 = @(x,y) exp(-((x-x0).^2+(y-y0).^2)/0.016);
fsrc = @(t,x,y) 0;
Tend = 1;
CFL = 0.5;
ht = Tend / 2048;
bct = 'dddd';
bcv = [0 0 0 0];
bp = [-1 1 -1 1];

[x,y,s,et]=adv_diff_2d([nx,ny],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct);

xi = zwgll(nx); yi = zwgll(ny);
xo = linspace(-1,1,64); yo = linspace(-1,1,64);
Jx = interp_mat(xo,xi); Jy = interp_mat(yo,yi);

[xx,yy] = ndgrid(xo, yo);
s_intp = tensor2(Jy, Jx, s);

figure(1, 'Units', 'inches', 'Position', [2 2 4 4]); box on;
    mesh(x,y,s,lw,1)
    view(-37.5, 20);
    zlim([0 0.4]);
    pbaspect([1 1 1]);
    xlabel('$x$',intp,ltx); ylabel('$y$',intp,ltx); zlabel('$u(x,y)$',intp,ltx);
    set(gca, lw, 1, fn, 'serif', fs, 12)
    savefig_pdf('02_01_adv_diff_mesh')
