addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib');
premeable;

% E2V = [1 2 3 4];
% VP = [0 0;
%      1 0;
%      0 1;
%      1 1;]
% read_quad_mesh(E2V, VP)
nex = 2; ney = 2;
[X,Y] = gen_box_2d(8, nex, ney);
% figure; hold on;
% axis equal
% for i=1:nex; for j=1:ney;
%     x = squeeze(X(:,(i-1)*ney+j,:));
%     y = squeeze(Y(:,(i-1)*ney+j,:));
%     mesh(x, y, x)
% end;end;
