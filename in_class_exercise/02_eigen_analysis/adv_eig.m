addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib')
premeable;

%% Eigen analysis for the advection problem

% 1D constant c
% c = 1;
% n = 100;
% [Ah,Bh,Ch,Dh,z,w] = semhat(n);
% R = speye(n+1)(2:end,:);    % Restriction
% B = R*Bh*R';
% Q = speye(n+1)(:,2:end); Q(1,n)=1;  % Periodic
% C = Q'*Ch*Q;
% [S,Lam] = eig(C,B);
% l = reshape(full(diag(Lam)), n, 1);
% j = sqrt(-1);
% k = -(n-1)/2:(n-1)/2;
% la = reshape(-1j * c * k * pi, n, 1);
% l = sort(imag(l));
% la = sort(imag(la));
% figure
% hold on
% % scatter(real(l), imag(l), 36)
% % scatter(real(la), imag(la), 36)
% plot(1:n,l,lw,1.5,1:n,la,lw,1.5)
% % semilogy(1:n,abs(l-la))
% hold off

% 2D plane rotation / straining flow with over integration
nxl = 19;
nyl = nxl;
nx = 18;
ny = nx;
[Ahx,Bhx,Chx,Dhx,zx,wx] = semhat(nx);
[Ahy,Bhy,Chy,Dhy,zy,wy] = semhat(ny);
[zxl, wxl] = zwgl(nxl);
[zyl, wyl] = zwgl(nyl);
Jx = interp_mat(zxl, zx);
Jy = interp_mat(zyl, zy);
Bxl = sparse(diag(wxl));
Byl = sparse(diag(wyl));


Ix = speye(nx+1); Iy = speye(ny+1);
Rx = Ix(2:end-1,:); Ry = Iy(2:end-1,:);
Bx = Rx*Bhx*Rx'; By = Ry*Bhy*Ry';
px = ones(length(zx),1); qx = zy;   % cx = y
py = -zx; qy = ones(length(zy),1);  % cy = -x
px = zx; qx = ones(length(zy),1);   % cx = x
py = ones(length(zx),1); qy = -zy;  % cy = -y
Px = sparse(diag(px)); Qx = sparse(diag(qx));
Py = sparse(diag(py)); Qy = sparse(diag(qy));

Cx = kron(Qx,Px); Cy = kron(Qy,Py);
R = kron(Ry,Rx);
B = sparse(kron(By, Bx));
Bl = kron(Byl, Bxl);
J = kron(Jy, Jx);
Clx = sparse(diag(J*diag(Cx)));
Cly = sparse(diag(J*diag(Cy)));

Cb = kron(Bhy*Qx,Bhx*Px*Dhx) + kron(Bhy*Qy*Dhy,Bhx*Py);
C = R*Cb*R';

Cbl = J'*Bl*(Clx*J*kron(Iy,Dhx) + Cly*J*kron(Dhy,Ix));
Cl = R*Cbl*R';

[V,D] = eig(C,B);
[Vl,Dl] = eig(Cl,B);
d = diag(D);
dl = diag(Dl);
scatter(real(dl), imag(dl), 50)