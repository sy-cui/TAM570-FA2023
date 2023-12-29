addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

function [R,T,E,it_t,res_t] = diffusion_convergence(N,ratio);
    N = [N,N,N];
    tol=1e-16; max_iter=1000;

    % Basic setup
    nx=N(1);ny=N(2);nz=N(3);
    [dim,x,w,D,wm,Jm,Jf] = set_mono_param(N,0,1.5,ceil(64/N(1)));
    [r,s,t] = x{:};
    x=r;y=s;z=t;
    [X,Y,Z] = ndgrid(x,y,z); 
    [X,Y,Z,Rad,Pol,Azi] = morph_sphere(X,Y,Z,ratio,1.0);
    Xf=t3w(Jf,X); Yf=t3w(Jf,Y); Zf=t3w(Jf,Z); 
    [Rx,G,B,Jac] = geom_elem_3D(X,Y,Z,D,w);
    vol = 4/3*pi*(1-ratio^3)
    [JRx,JD,Bm] = set_dealiase_op(dim,D,wm,Jm,Rx,Jac);

    % Time steps
    dl = {2,2,2};

    Rt=set_restriction(N,{'p','p','n','n','d','d'}); [At,St,Lt]=neumann_op(dim,Rt,w,D,dl); 

    % Solution fields
    T = zeros(N+1);
    T(:,:,1) = 1; Tb = T;
    Lt_inv = 1.0 ./ (Lt);
    Pinv_t = 1;

    Th = -viscous_op(Tb,{1,1,1},D,G,B,0,1);
    [T,it_t,res_t,~] = pcg(t3w(Rt,Th),Pinv_t,St,Lt_inv,Rt,D,G,B,0,1,max_iter,tol);
    T = T + Tb;

    R = squeeze(t3w(Jf,Rad)(1,1,:));
    T = squeeze(t3w(Jf,T)(1,1,:));
    E = ratio / (1 - ratio) * (1./ R - 1);

end;

[R,T,E,it,res] = diffusion_convergence(16,0.5);

figure(1, 'Units', 'inches', 'Position', [2 2 6 6])
    plot(R,T,'o-k',lw,1.5,'DisplayName','Simulation')
    hold on;
    plot(R,E,'-r',lw,1.5,'DisplayName','Exact solution')
    hold off;
    xlabel('Radial distance $r$', intp, ltx);
    ylabel('Temperature $T$', intp, ltx);
    legend(intp, ltx);
    set(gca, fs, 16, fn, 'serif')
    savefig_png('figures/radial_temperature')
    

