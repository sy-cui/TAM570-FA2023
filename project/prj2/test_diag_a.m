addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% 2D
N = [5,10];
Dhx=rand(N(1)); Dhy=rand(N(2));
Ix = eye(N(1)); Iy = eye(N(2));
G_cell = cell(1,3);
for i = 1:3; G_cell{i} = rand(N); end;
dA = diag_a_2d(Dhx,Dhy,G_cell{1},G_cell{2},G_cell{3});

Grr = diag(reshape(G_cell{1},[],1));
Gss = diag(reshape(G_cell{2},[],1));
Grs = diag(reshape(G_cell{3},[],1));
Dx=kron(Iy, Dhx);
Dy=kron(Dhy, Ix);
D = [Dx; Dy];
G = [Grr Grs; Grs Gss];
A = D' * G * D;
dA_e = reshape(diag(A), N);
sum(sum(abs(dA - dA_e)))

% 3D
N = [6,7,8];
Dhx=rand(N(1)); Dhy=rand(N(2)); Dhz=rand(N(3));
Ix = eye(N(1)); Iy = eye(N(2)); Iz = eye(N(3));
G_cell = cell(1,6);
for i = 1:6; G_cell{i} = rand(N); end;
dA = diag_a_3d({Dhx,Dhy,Dhz},G_cell);

Grr = diag(reshape(G_cell{1},[],1));
Gss = diag(reshape(G_cell{2},[],1));
Gtt = diag(reshape(G_cell{3},[],1));
Grs = diag(reshape(G_cell{4},[],1));
Grt = diag(reshape(G_cell{5},[],1));
Gst = diag(reshape(G_cell{6},[],1));
Dx=kron3(Iz,Iy,Dhx);
Dy=kron3(Iz,Dhy,Ix);
Dz=kron3(Dhz,Iy,Ix);
D = [Dx; Dy; Dz];
G = [Grr Grs Grt; Grs Gss Gst; Grt Gst Gtt];
A = D' * G * D;
dA_e = reshape(diag(A), N);
sum(sum(sum(abs(dA - dA_e))))