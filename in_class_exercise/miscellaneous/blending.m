addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;


Xd = [
    0.5 0.7 0.6 0.4;
    1.2 1.4 1.3 1.0; 
    2.4 2.4 2.3 2.2;
    3.0 2.9 2.8 2.5
];
Yd = [
    0.5 1.0 1.7 2.1;
    0.6 1.1 1.8 2.4;
    0.8 1.3 2.0 2.4;
    1.1 1.5 2.0 2.3;
];

[X, Y] = laplace_blend_2d(Xd, Yd);

[nx, ny] = size(Xd);
[xi, ~] = zwgll(nx-1); [yi, ~] = zwgll(ny-1); 
[xo, ~] = zwuni(32); [yo, ~] = zwuni(32);
Jx = interp_mat(xo, xi); Jy = interp_mat(yo, yi);

Xf = Jx * X * Jy'; Yf = Jx * Y * Jy';

figure; hold on
scatter(Xd, Yd, 'k', 'o')
plot(Xf, Yf, '-r', Xf', Yf', '-r')
hold off



