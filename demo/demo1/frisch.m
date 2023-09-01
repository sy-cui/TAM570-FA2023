format compact

%% Frisch Map.  Usage:  n=60; frisch

u=zeros(n,1); r=zeros(n,1);

v=0.11110101010101; w=0.11110101010102; d0=abs(v-w);

t=1:n;
for j=1:n;
  v=1-2*v*v; u(j)=v;
  w=1-2*w*w; r(j)=w;
  s=[ j v w w-v ];
  z=sprintf('%4d %16.11f  %16.11f  %16.12f',s);
  disp(z)
end;


plot(t,u,'ro-',t,r,'bo-'); axis square
xlabel('Iteration Number, j','fontsize',16);
ylabel('Solution: v_j and w_j','fontsize',16);
title('Evolution of Frisch Map, different ICs','fontsize',16);
set (gcf,'color',[1 1 1]) % set the background color to white
set(gcf,'PaperUnits','normalized');set(gcf,'PaperPosition',[0 0 1 1])
print -dpng 'frisch_soln.png'
% !open frisch_soln.png

model = d0*2.^t; model=min(2,model);
semilogy(t,abs(r-u),'g.-',t,model,'linewidth',1.2); axis square
xlabel('Iteration Number, j','fontsize',16);
ylabel('Difference: |v_j-w_j|','fontsize',16);
title('Separation d_j=|v_j-w_j| and  2^jd_0 vs. j','fontsize',16);
set (gcf,'color',[1 1 1]) % set the background color to white
set(gcf,'PaperUnits','normalized');set(gcf,'PaperPosition',[0 0 1 1])
print -dpng 'frisch_growth.png'
% !open frisch_growth.png
