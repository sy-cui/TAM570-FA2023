addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

H = 1; alpha = 1.0; L = alpha*2*pi*H; epsilon = 1e-5;

nv = [32 32]; mv = ceil(nv * 1.5);
Re = 7500; w_i = 0.002234975649;
Pr = 5e-4;
Tend = 200; CFL=0.5; Vscale=1;
ht = 0;
bp = [0 L -H H];

u_ini = @(x,y) (1-y.^2)+epsilon*cos(2*pi*x/L).*sin(pi*y/H);
v_ini = @(x,y) 0*x;
t_ini = @(x,y) 0*x;

u_src = @(t,x,y) 2/Re+0*x;
v_src = @(t,x,y) 0*x;
t_src = @(t,x,y) 0*x;

u_ext = @(t,x,y) (1-y.^2);
v_ext = @(t,x,y) 0*x;
t_ext = @(t,x,y) 0*x;

f0 = {u_ini, v_ini, t_ini};
fq = {u_src, v_src, t_src};
fe = {true, u_ext, true, v_ext, false, t_ext};

[x,y,soln,err,et] = ns2d(nv,mv,Re,Pr,Tend,CFL,Vscale,ht,bp,f0,fq,fe);
t = linspace(0.0, Tend, length(err{4})+1)(2:end);
perturb = err{4};
save("orr_sommerfield.mat", "t", "perturb")

load("orr_sommerfield.mat")
perturb = perturb * sqrt(4*pi);
figure(2, 'Units', 'inches', 'Position', [2 2 6 6]); box on;
hold on;
coeff = 6.5e-6;
semilogy(t, perturb, '-k', lw, 1.5);
semilogy(t, coeff*exp(w_i*t), '-r', lw, 1.5);
hold off;
xlabel("Time");
ylabel("$L_2$ Perturbation magnitude", intp, ltx);
legend("Numerical difference", "Model $\\alpha e^{\\omega_i t}$", intp, ltx)
title(sprintf("$\\alpha=%d, N_x=%d, N_y=%d$", coeff, nv(1), nv(2)), intp, ltx)
set(gca, fs, 16, fn, "serif", lw, 1.5)
savefig_pdf("orr_sommerfield_perturb")


