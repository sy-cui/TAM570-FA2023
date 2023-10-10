%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

run_spatial_convergence = false;
run_temporal_convergence = false;

nu = 0.01;
fcx = @(x,y) 0*x;
fcy = @(x,y) 0*x;
mode = 1.5;
f0 = @(x,y) cos(mode*pi*x).*cos(mode*pi*y);

% Analytical solution
tau = 0.5 / nu / (mode * pi)^2;
fsoln = @(t,x,y) exp(-t/tau)*f0(x,y);

fsrc = @(t,x,y) 0;
Tend = log(2) * tau;
CFL = 0.5;
bct = 'dddd';
bcv = [0 0 0 0];
bp = [-1 1 -1 1];
fa = {true, fsoln};

if run_spatial_convergence;
    % Spatial convergence
    ns = [4:2:32 36:4:48];
    espc = ns * 0;
    for i=1:length(ns);
        N = ns(i);
        fprintf("Running N = %d \n", N);
        [x,y,s,et]=adv_diff_2d([N,N],nu,fcx,fcy,f0,fsrc,Tend,CFL,Tend/8192,bct,bcv,bp,fa);
        soln = fsoln(Tend,x,y);
        espc(i) = max(max(soln - s));
    end;
    save('diffusion_spatial_convergence.mat', 'ns', 'espc')
end;

if run_temporal_convergence;
    % Temporal convergence
    tsteps = 2.^(2:13);
    etpr = tsteps * 0;
    for i=1:length(tsteps);
        ht = Tend / tsteps(i);
        fprintf("Running tsteps = %d \n", tsteps(i));
        [x,y,s,et]=adv_diff_2d([128,128],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct,bcv,bp,fa);
        soln = fsoln(Tend,x,y);
        etpr(i) = max(max(soln - s));
    end;
    save('diffusion_temporal_convergence.mat', 'tsteps', 'etpr')
end;

load('diffusion_spatial_convergence.mat')
load('diffusion_temporal_convergence.mat')
figure(1, 'Units', 'inches', 'Position', [2 2 9 4]);
    subplot(1, 2, 1); box on;
    semilogy(ns, espc, '-ok', lw, 1.5);
    xlabel('$N$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    set(gca, lw, 1, fs, 12, fn, 'serif')

    subplot(1, 2, 2); box on;
    hold on;
    loglog(Tend./tsteps, etpr, '-ok', lw, 1.5);
    loglog(Tend./tsteps(1:end-1), 5e-2*(Tend./tsteps)(1:end-1).^3, '--r', lw, 1.5)
    hold off;
    xlabel('$\Delta t$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    set(gca, lw, 1, fs, 12, fn, 'serif')
    savefig_pdf('02_02_convergence')





