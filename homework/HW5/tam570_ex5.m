addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
nx = 40;
ny = nx;
T = 2*pi;

[Ahx,Bhx,Chx,Dhx,zx,wx] = semhat(nx);
[Ahy,Bhy,Chy,Dhy,zy,wy] = semhat(ny);
[R,S] = ndgrid(zx,zy);

dt = 0.5*(zx(2)-zx(1))/sqrt(2);
nt = ceil(T/dt);
dt = T/nt;

u0 = ones(nx+1,ny+1);          % initial condition
x0 = 0.5; y0 = 0.0;

Ihx = speye(nx+1); Rx=speye(nx+1); Rx=Rx(2:end-1,:);
Ihy = speye(ny+1); Ry=speye(ny+1); Ry=Ry(2:end-1,:);
velx = zeros(nx+1,ny+1);
vely = zeros(nx+1,ny+1);

for i = 1:nx+1
    for j = 1:ny+1
        u0(i,j) = exp(-((zx(i) - x0)^2 + (zy(j) - y0)^2) * (1 / 0.016));
        velx(i,j) = -zy(j);
        vely(i,j) = zx(i);
    end
end

Ax = Rx*Ahx*Rx'; Ay = Ry*Ahy*Ry';
Bx = Rx*Bhx*Rx'; By = Ry*Bhy*Ry';
Bix  = diag(1 ./ diag(Bx)); Biy = diag(1 ./ diag(By));
Cx = Rx*Chx*Rx'; Cy = Ry*Chy*Ry';
%Dx = Rx*Dhx*Rx'; Dy = Ry*Dhy*Ry';
Ix = Rx*Ihx*Rx'; Iy = Ry*Ihy*Ry';
velx = Rx * velx * Rx'; vely = Ry * vely * Ry';

u_arr = zeros((nx+1),(ny+1),nt+2);
u_arr(:,:,1) = u0;
u_arr(:,:,2) = u0;
u_arr(:,:,3) = u0;

[Sx, Lamx] = gen_eig_decomp(Cx,Bx);
[Sy, Lamy] = gen_eig_decomp(Cy,By);

%lamx = reshape(Lamx, [length(Lamx), 1])
%lx = full(diag(Lamx));
%Lam = 1 ./ reshape((reshape(Lamx, length(Lamy), []) + Lamy), (ny+1)*(nx+1), (ny+1)*(nx+1));

lx = full(diag(Lamx)); ly = full(diag(Lamy));

p1 = reshape(lx, [length(lx), 1]);       
p2 = reshape(ly, [1, length(ly)]);
Lam = 1 ./ (p1 + p2);

bk = [1.5 -2 0.5];
ak = [2 -1];

bk = [11 -18 9 -2] / 6;
ak = [3 -3 1];

for i=4:(nt+2)
    f1 = tensor2(Ry, Rx, u_arr(:,:,i-1));
    f2 = tensor2(Ry, Rx, u_arr(:,:,i-2));
    f3 = tensor2(Ry, Rx, u_arr(:,:,i-3));
    ext = (ak(1) .* f1 + ak(2) .* f2 + ak(3) .* f3);
    Vx = velx .* tensor2(By, Cx, ext);
    Vy = vely .* tensor2(Cy, Bx, ext);
    rhs = - bk(2) .* f1 - bk(3) .* f2 - bk(4) .* f3 - tensor2(Biy, Bix, dt .* (Vx + Vy));
    uh = rhs / bk(1);
    u_arr(:,:,i) = tensor2(Ry', Rx', uh);

    if mod(i,100) == 0;
        figure(1)
        surf(R, S, u_arr(:,:,i));
        pause(0.01)
    end;

    %hl = ( I + 0.5 .* dt .* alpha .*L);
    %hr = ( I - 0.5 .* dt .* alpha .*L);
    %u_arr(:,i) = R'*(inv(hl) * (hr) * (R * u_arr(:,i-1)));
end
% surf(R, S, u_arr(:,:,10));