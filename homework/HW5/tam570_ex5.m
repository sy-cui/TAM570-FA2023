addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
narr = 2:80;
cfl = [0.25 0.5 1];
err = zeros(length(cfl), length(narr));

for l=1:length(cfl)
    for k=1:length(narr)
        nx = narr(k);
        ny = nx;
        T = 2 * pi;
        nu = 0.00;
        CFL = cfl(l);
        
        [Ahx,Bhx,Chx,Dhx,zx,wx] = semhat(nx);
        [Ahy,Bhy,Chy,Dhy,zy,wy] = semhat(ny);
        [R,S] = ndgrid(zx,zy);
        
        dt = CFL*(zx(2)-zx(1))/1;
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
                %u0(i,j) = cos(pi * zx(i) * 0.5) * cos(pi * zy(j) * 0.5);
                velx(i,j) = 1 * -zy(j);
                vely(i,j) = 1 * zx(i);
            end
        end
        
        Ax = Rx*Ahx*Rx'; Ay = Ry*Ahy*Ry';
        Bx = Rx*Bhx*Rx'; By = Ry*Bhy*Ry';
        Bix  = diag(1 ./ diag(Bx)); Biy = diag(1 ./ diag(By));
        Cx = Rx*Chx*Rx'; Cy = Ry*Chy*Ry';
        %Dx = Rx*Dhx*Rx'; Dy = Ry*Dhy*Ry';
        Ix = Rx*Ihx*Rx'; Iy = Ry*Ihy*Ry';
        velx = Rx * velx * Rx'; vely = Ry * vely * Ry';
        
        u_arr = zeros((nx+1),(ny+1),nt+3);
        u_arr(:,:,1) = u0;
        u_arr(:,:,2) = u0;
        u_arr(:,:,3) = u0;
        
        [Sx, Lamx] = gen_eig_decomp(Ax, Bx);
        [Sy, Lamy] = gen_eig_decomp(Ay, By);
        
        %lamx = reshape(Lamx, [length(Lamx), 1])
        %lx = full(diag(Lamx));
        %Lam = 1 ./ reshape((reshape(Lamx, length(Lamy), []) + Lamy), (ny+1)*(nx+1), (ny+1)*(nx+1));
        
        lx = full(diag(Lamx)); ly = full(diag(Lamy));
        
        bk = [11 -18 9 -2] / 6;
        ak = [3 -3 1];
        
        p1 = reshape(lx, [length(lx), 1]);       
        p2 = reshape(ly, [1, length(ly)]);
        Lam = 1 ./ ((p1 + p2) .* (dt) .* nu + bk(1));
        
        for i=4:nt+3
            f1 = tensor2(Ry, Rx, u_arr(:,:,i-1));
            f2 = tensor2(Ry, Rx, u_arr(:,:,i-2));
            f3 = tensor2(Ry, Rx, u_arr(:,:,i-3));
            ext = (ak(1) .* f1 + ak(2) .* f2 + ak(3) .* f3);
            Vx = velx .* tensor2(By, Cx, ext);
            Vy = vely .* tensor2(Cy, Bx, ext);
            rhs = - bk(2) .* tensor2(By, Bx, f1) - bk(3) .* tensor2(By, Bx, f2) - bk(4) .* tensor2(By, Bx, f3) - dt .* (Vx + Vy);
            un = tensor2(Sy, Sx, Lam .* tensor2(Sy', Sx', rhs));
            %un = (1 / (bk(1))) .* un;
            u_arr(:,:,i) = tensor2(Ry', Rx', un);
            %u_arr(:,:,1:3) = u_arr(:,:,2:4);
        end
        uan = u0;
        %tau = 2 / ( pi * pi * nu);
        %uan = exp(-T/tau) * u0;
        err(l,k) = max(abs(u_arr(:,:,end) - uan),[],'all');
        figure(2);
        surf(R, S, u_arr(:,:,end));
        xlabel("$X$", Interpreter="latex");
        ylabel("$Y$", Interpreter="latex");
        zlabel('$u(x,y,t)|_{t = T}$', 'Interpreter','latex')
        pause(0.01)
        hAxes.TickLabelInterpreter = 'latex';
        title(['N = ', num2str(narr(k) + 1)], Interpreter="latex")
        pause(0.01)
    end
end
