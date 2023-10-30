addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

run_case = false;
plot_nusselt = false;
plot_temp = true;

% Benchmark comparison with Ouertatani et al. (2008)
nv = [32, 32];
Ra = 1e5;
Pr = 0.71;
Tend = 2;
CFL = 0.005;
epsilon = 0;
kc = 2 * pi;
x_bc = "d";
uy_bc = "DD";

[xu, ~] = zwuni(128); xu = 0.5 * (xu + 1); 

if run_case;
    [x,y,U,V,T,Nu,ht,et,~] = rayleigh_benard(nv,Ra,Pr,Tend,CFL,epsilon,kc,x_bc,uy_bc,1e-5);
    save("data/Ra5_alt.mat", "x", "y", "U", "V", "T", "ht", "Nu", "et")
end; % run_case

% data = dlmread('nusselt.csv');
% Ra4_x = data(:, 1); Ra4_y = data(:, 2); 
% Ra5_x = data(:, 3); Ra5_y = data(:, 4);
% Ra6_x = data(:, 5); Ra6_y = data(:, 6);

% save('data/nusselt_data.mat', 'Ra4_x', 'Ra4_y', 'Ra5_x', 'Ra5_y', 'Ra6_x', 'Ra6_y')
% load ('data/nusselt_data.mat')

if plot_nusselt;
    figure(2, 'Units', 'inches', 'Position', [2 2 10 6]); box on;
    hold on
    % scatter(Ra4_x, Ra4_y, 'k', 'DisplayName', '$Ra=10^4$')
    % scatter(Ra5_x, Ra5_y, 'r', 'DisplayName', '$Ra=10^5$')
    % scatter(Ra6_x, Ra6_y, 'b', 'DisplayName', '$Ra=10^6$')

    load('data/Ra3.mat')
    Jx = interp_mat(xu, x); Nu = Jx * Nu;
    [maxNu, idx] = max(Nu);
    fprintf("Maximum Nu for Ra=1e3 is %d occurring at %d\n", maxNu, xu(idx));
    plot(xu, Nu, '-k', lw, 2, 'DisplayName', '$Ra=10^3$')

    load('data/Ra4.mat')
    Jx = interp_mat(xu, x); Nu = Jx * Nu;
    [maxNu, idx] = max(Nu);
    fprintf("Maximum Nu for Ra=1e4 is %d occurring at %d\n", maxNu, xu(idx));
    plot(xu, Nu, '-r', lw, 2, 'DisplayName', '$Ra=10^4$')

    load('data/Ra5.mat')
    Jx = interp_mat(xu, x); Nu = Jx * Nu;
    [maxNu, idx] = max(Nu);
    fprintf("Maximum Nu for Ra=1e5 is %d occurring at %d\n", maxNu, xu(idx));
    plot(xu, Nu, '-b', lw, 2, 'DisplayName', '$Ra=10^5$')

    load('data/Ra6.mat')
    Jx = interp_mat(xu, x); Nu = Jx * Nu;
    [maxNu, idx] = max(Nu);
    fprintf("Maximum Nu for Ra=1e6 is %d occurring at %d\n", maxNu, xu(idx));
    plot(xu, Nu, '-g', lw, 2, 'DisplayName', '$Ra=10^6$')
    hold off

    xlim([0 1]); ylim([0 12]);
    xlabel("$x$", intp, ltx); ylabel("$Nu_L$", intp, ltx);
    set(gca, fs, 32, fn, 'serif', lw, 1.5)
    legend(intp, ltx, fs, 26, 'location', 'northwest');

    savefig_png('figures/nusselt_profile')

end; % plot_nusselt

%% Temperature profile
if plot_temp;
    figure(2, 'Units', 'inches', 'Position', [5 5 8 8])
    box on; hold on;
        load('data/Ra3.mat')
        [~, w] = zwgll(length(x) - 1); B = diag(0.5 * w); 
        J = interp_mat(xu, x); Tave = J * sum(B * T, 1)';
        plot(xu, Tave, '-k', lw, 2, 'DisplayName', '$Ra=10^3$')

        load('data/Ra4_alt.mat')
        [~, w] = zwgll(length(x) - 1); B = diag(0.5 * w); 
        J = interp_mat(xu, x); Tave = J * sum(T * B, 2);
        plot(xu, Tave, '-r', lw, 2, 'DisplayName', '$Ra=10^4$')

        load('data/Ra5_alt.mat')
        [~, w] = zwgll(length(x) - 1); B = diag(0.5 * w); 
        J = interp_mat(xu, x); Tave = J * sum(T * B, 2);
        plot(xu, Tave, '-b', lw, 2, 'DisplayName', '$Ra=10^5$')

        load('data/Ra6_alt.mat')
        [~, w] = zwgll(length(x) - 1); B = diag(0.5 * w); 
        J = interp_mat(xu, x); Tave = J * sum(T * B, 2);
        plot(xu, Tave, '-g', lw, 2, 'DisplayName', '$Ra=10^6$')
    hold off

        xlim([0 1]); ylim([0 1]);
        xlabel('$x/L_x$', intp, ltx); ylabel('$\bar{T}$', intp, ltx);
        set(gca, fs, 32, fn, 'serif', lw, 1.5)
        legend(intp, ltx, fs, 20);
        savefig_png('figures/temp_profile')

end; % plot_temp
