addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

run = false;

nv = [32, 32];
Pr = 5;
Tend = 300;
CFL = 0.5;
epsilon = 1e-3;
x_bc = "p";

ht = 3.4718e-03;

Ra_arr = [665:5:695]; uy_bc = "NN"; kc = 2.2214;
% Ra_arr = [1110:5:1140]; uy_bc = "DN"; kc = 2.682;
% Ra_arr = [1720:10:1780]; uy_bc = "DD"; kc = 3.117;

if run;
    for Ra = Ra_arr;
        fprintf("Running Ra = %d\n", Ra);
        fname = strcat('data/Ra', num2str(Ra), '_', uy_bc, '.mat');
        [x,y,U,V,T,Nu,ht,et,KE] = rayleigh_benard(nv,Ra,Pr,Tend,CFL,epsilon,kc,x_bc,uy_bc);
        save(fname,'x','y','U','V','T','Nu','ht','et','KE');
    end;
end;

figure(1, 'Units', 'inches', 'Position', [2 2 6 6]);
box on; hold on;

KE_arr = Ra_arr * 0;
GR = KE_arr + 0;
for k = 1:length(Ra_arr);
    Ra = Ra_arr(k);
    fname = strcat('data/Ra', num2str(Ra), '_', uy_bc, '.mat');
    load(fname); KE_arr(k) = KE(end);

    label = sprintf('$Ra = %d$', Ra);
    semilogy([1:length(KE)]*ht, KE, lw, 2, 'DisplayName', label)
    % cp = ceil(length(KE)/2);
    % LinKE = log(KE(cp - 100 : cp + 100));
    % t = [1:length(LinKE)] * ht;
    % growth_rate = polyfit(t, LinKE, 1); 
    % GR(k) = growth_rate(1);
end;
hold off
% legend(intp, ltx, 'location', 'southeast')
% xlabel('Time'); ylabel('$E_k$', intp, ltx);
% set(gca, fs, 24, lw, 2, fn, 'serif')
% savefig_png('figures/kinetic_energy_evolution')


figure(2, 'Units', 'inches', 'Position', [2 2 6 6]);
box on; hold on;
    p = polyfit(Ra_arr, KE_arr, 1);
    Ra_c = -p(2) / p(1);
    Ra_c
    Ra_arr
    KE_arr
    pause
    x = linspace(Ra_c, Ra_arr(end) + 5, 10);
    y = p(1) * x + p(2);

    plot(x, y, '--k', lw, 1.5);
    scatter(Ra_arr, KE_arr, 80, 'k', 'd', lw, 1.5);
    scatter([Ra_c], [0], 80, 'r', '*', lw, 1.5);

    xlabel('$Ra$', intp, ltx); ylabel("$\\bar{E}_k$", intp, ltx);
    set(gca, fs, 24, lw, 2, fn, 'serif');
    savefig_png("figures/ke_linfit")

hold off;
