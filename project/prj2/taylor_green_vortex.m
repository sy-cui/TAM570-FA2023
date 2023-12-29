addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

Utrans = 1; Vtrans = 1; wn = pi;

nv = [16 16]; mv = [24 24];
Re = 100; Pr = 5e-4;
Tend = 1; CFL=0.1; Vscale=1;
ht = 0;
bp = [-1 1 -1 1];

u_ini = @(x,y) sin(wn*x).*cos(wn*y) + Utrans;
v_ini = @(x,y) -cos(wn*x).*sin(wn*y) + Vtrans;
t_ini = @(x,y) 0*x;

u_src = @(t,x,y) 0*x;
v_src = @(t,x,y) 0*x;
t_src = @(t,x,y) 0*x;

u_ext = @(t,x,y) sin(wn*(x-Utrans*t)).*cos(wn*(y-Vtrans*t))*exp(-2*wn*wn*t/Re) + Utrans;
v_ext = @(t,x,y) -cos(wn*(x-Utrans*t)).*sin(wn*(y-Vtrans*t))*exp(-2*wn*wn*t/Re) + Vtrans;
t_ext = @(t,x,y) 0*x;

f0 = {u_ini, v_ini, t_ini};
fq = {u_src, v_src, t_src};
fe = {true, u_ext, true, v_ext, false, t_ext};

ns = [32];
er = ns * 0;
for i = 1:length(ns);
    N = ns(i);
    M = ceil(1.5*N);
    [x,y,soln,err,et] = ns2d([N,N],[M,M],Re,Pr,Tend,CFL,Vscale,ht,bp,f0,fq,fe);
    er(i) = err{4}(end);
end;

% semilogy(ns, er, '-k', lw, 1.5);
