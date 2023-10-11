%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

run_convergence = false;

x0 = 0.5; y0 = 0;
nu = 0;
fcx = @(x,y) -y;
fcy = @(x,y) x;
f0 = @(x,y) exp(-((x-x0).^2+(y-y0).^2)/0.016);
fsrc = @(t,x,y) 0;
Tend = 2*pi;
CFL = 0.5;
ht = 0;
bct = 'dddd';
bcv = [0 0 0 0];
bp = [-1 1 -1 1];
fa = {true, @(t,x,y) exp(-((x-x0*cos(t)).^2+(y-x0*sin(t)).^2)/0.016)};
% [x,y,s,et]=adv_diff_2d([40,40],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct);
% fsoln = @(x,y) exp(-((x-x0*cos(Tend)).^2+(y-x0*sin(Tend)).^2)/0.016);

if run_convergence;
    % Convergence
    n = [4:2:60 64:4:80];
    errors = n*0;
    for i=1:length(n);
        N = n(i);
        fprintf('Solving n = %d \t', N);
        [x,y,s,et]=adv_diff_2d([N,N],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct,bcv,bp,fa);
        errors(i) = max(max(abs(s-f0(x,y))));
        fprintf('Elapsed time is %f \n', et);
    end;

    save(strcat(['adv_errs_cfl_',strrep(num2str(CFL),'.','pt'),'.mat']),'CFL','n','errors')
end;

figure(1, 'Units', 'inches', 'Position', [2 2 8 5]); box on;
    hold on
    load('adv_errs_cfl_0pt25.mat')
    semilogy(n, errors, '-ok', lw, 1.5)
    load('adv_errs_cfl_0pt5.mat')
    semilogy(n, errors, '-or', lw, 1.5)
    load('adv_errs_cfl_1.mat')
    semilogy(n, errors, '-ob', lw, 1.5)
    xlabel('$N$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    plot(n, 80*exp(-n/3.5), '-.', lw, 1.5, 'color', [0 0.4470 0.7410]);
    plot(n, 7e5*n.^(-6), '-.', lw, 1.5, 'color', [0.8500 0.3250 0.0980]);
    plot(n, 7/8*1e5*n.^(-6), '-.', lw, 1.5, 'color', [0.4940 0.1840 0.5560]);
    ylim([1e-7 10])
    legend('$\text{CFL}=0.25$', '$\text{CFL}=0.5$', '$\text{CFL}=1.0$', 
        '$E = e^{-\mu N}$', '$E = \alpha N^{-6}$', '$E = \alpha N^{-6}/8$', intp, ltx);
    hold off
    set(gca,fs,16,fn,'serif',lw,1.5);
    savefig_pdf('01_01_convergence')

% Mesh plots
CFL = 0.5;
figure(2, 'Units', 'inches', 'Position', [2 2 12 4]);
n = [10 30 60];
for i=1:3;
    subplot(1,3,i); box on;
    N = n(i);
    [x,y,s,et]=adv_diff_2d([N,N],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct);
    mesh(x,y,s,lw,1);
    view(-37.5, 20);
    zlim([-0.2 1.0])
    xlabel('$x$',intp,ltx); ylabel('$y$',intp,ltx); zlabel('$u(x,y)$',intp,ltx);
    title(sprintf('$N = %d$',N), intp, ltx);
    set(gca, lw, 1, fn, 'serif', fs, 12)
end;
savefig_pdf('01_02_mesh')
