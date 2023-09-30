%% Setup
addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% Eigen analysis of Graetz problem via Legendre approximation
ny = 256;
alpha = 0.01;

[Ah, Bh, ~, ~, y, ~] = semhat(ny);
R = speye(ny + 1)(2:end-1, :);
ux = 1 - y.^2;

A = alpha*R*Ah*R';
B = R*(Bh.*ux)*R'; 
[S, Lam] = gen_eig_decomp(A, B);

% Comparison
k = 1:4;
S_diff = sin(0.5*pi*(y+1)*k);

figure(1, 'Units', 'inches', 'Position', [2 2 10 5])
subplot(1, 2, 1); hold on; box on;
subplot(1, 2, 2); hold on; box on;
for i=k;
    subplot(1, 2, 1);
    plot(y, S_diff(:,i), lw, 1)
    subplot(1, 2, 2);
    plot(y, [0; -S(:, i); 0], lw, 1)
end;
subplot(1, 2, 1); hold off;
subplot(1, 2, 2); hold off;
