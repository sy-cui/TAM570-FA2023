%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

x0 = 0.5; y0 = 0;
nu = 0;
fcx = @(x,y) -y;
fcy = @(x,y) x;
f0 = @(x,y) exp(-((x-x0).^2+(y-y0).^2)/0.016)+1;
fs = @(t,x,y) 0;
Tend = 2*pi;
CFL = 0.5;
ht = 0;
bct = 'dddd';
bcv = [1 1 1 1];

[x,y,s,et]=adv_diff_2d([80,80],nu,fcx,fcy,f0,fs,Tend,CFL,ht,bct,bcv);

% Convergence
% n = [4:2:60 64:4:80];
% errors = n*0;
% for i=1:length(n);
%     N = n(i);
%     fprintf('Solving n = %d \t', N);
%     [x,y,s,et]=adv_diff_2d([N,N],nu,fcx,fcy,f0,fs,Tend,CFL,ht,bct);
%     errors(i) = max(max(abs(s-f0(x,y))));
%     fprintf('Elapsed time is %f \n', et);
% end;

% semilogy(n, errors)
