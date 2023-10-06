%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

x0 = 0.5; y0 = 0;
nx = 80;
ny = 80;
nu = 0.1;
fcx = @(x,y) -cos(0.5*pi*x).*sin(0.5*pi*y);
fcy = @(x,y) sin(0.5*pi*x).*cos(0.5*pi*y);
% f0 = @(x,y) exp(-((x-x0).^2+(y-y0).^2)/0.016);
f0 = @(x,y) cos(0.5*pi*x).*cos(0.5*pi*y)

% Analytical solution
tau = 2 / (nu * pi^2);
fa = @(t,x,y) exp(-t/tau)*f0(x,y);

fs = @(t,x,y) 0;
Tend = log(2) * tau;
CFL = 0.5;
ht = 0;
bct = 'dddd';

[x,y,s,et]=adv_diff_2d([nx,ny],nu,fcx,fcy,f0,fs,Tend,CFL,ht,bct);
soln = fa(Tend,x,y);
max(max(abs(soln-s)))

