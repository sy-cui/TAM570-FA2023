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
[x,y,s,et]=adv_diff_2d([40,40],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct);

if run_convergence;
    % Convergence
    n = [4:2:60 64:4:80];
    errors = n*0;
    for i=1:length(n);
        N = n(i);
        fprintf('Solving n = %d \t', N);
        [x,y,s,et]=adv_diff_2d([N,N],nu,fcx,fcy,f0,fsrc,Tend,CFL,ht,bct);
        errors(i) = max(max(abs(s-f0(x,y))));
        fprintf('Elapsed time is %f \n', et);
    end;

    save(strcat(['adv_errs_cfl_',strrep(num2str(CFL),'.','pt'),'.mat']),'CFL','n','errors')
end;

figure(1, 'Units', 'inches', 'Position', [2 2 4 4]); box on;
    hold on
    load('adv_errs_cfl_0pt25.mat')
    semilogy(n, errors, '-ok', lw, 1.5)
    load('adv_errs_cfl_0pt5.mat')
    semilogy(n, errors, '-or', lw, 1.5)
    load('adv_errs_cfl_1.mat')
    semilogy(n, errors, '-ob', lw, 1.5)
    xlabel('$N$', intp, ltx); ylabel('$L_\infty$ error', intp, ltx);
    ylim([1e-7 10])
    legend('$CFL=0.25$', '$CFL=0.5$', '$CFL=1.0$', intp, ltx);
    hold off
    set(gca, fs, 12, fn, 'serif', lw, 1.5);
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
