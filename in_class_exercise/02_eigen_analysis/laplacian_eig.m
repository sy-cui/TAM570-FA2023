addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib')
premeable;

%% Laplacian operator eigenanalysis. 
k = 1:60;
L = 1;
lam_soln = k.^2 * pi^2 / L^2;

% Finite difference
N = k(end);
dx = L / N;
lam_fd = 2 / dx^2 * (1 - cos(pi*dx*k/L));

% Legendre spectral
N = k(end) + 1;
[Ah,Bh,Ch,Dh,z,w] =  semhat(N);
Ih=eye(N+1);

R=Ih(2:N,:);

A=R*Ah*R'*(2/L);
B=R*Bh*R'*(L/2);

[S,Lam]=gen_eig_decomp(A,B);
lam_gll = diag(full(Lam));
hold on
semilogy(k, abs(lam_soln(1:end)-lam_gll'), '-k', lw, 1)
xline(2/pi*N, '--r')
hold off