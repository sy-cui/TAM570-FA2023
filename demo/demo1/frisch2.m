lw='linewidth';               %%% Plotting defs
fs='fontsize';                %%% Plotting defs
intp = 'interpreter';         %%% Plotting defs
ltx  = 'latex';               %%% Plotting defs
format compact;

close all


%% Frisch Map.  Usage:  n=60; frisch

[zf,wf]=zwuni(1000);

n=60; u=zeros(n+1,1); r=zeros(n+1,1);

v=0.11110101010101; w=0.11110101010102; d0=abs(v-w);
u(1)=v;
r(1)=w;

t=1:n+1;
for j=1:n;
  vo=v; v=1-2*v*v; u(j+1)=v;
  wo=w; w=1-2*w*w; r(j+1)=w;
  s=[ j v w w-v ];
  z=sprintf('%4d %16.11f  %16.11f  %16.12f',s);
  disp(z)

  hold off; plot(zf,0*zf,'k-',zf,1-2*zf.*zf,'g-',lw,2); axis square; hold on;
  xlabel('x-value input',fs,14); ylabel('y-value output',fs,14);
  title('Frisch Nonlinear Map: v_n = 1-v_{n-1}^2',fs,14);
  plot(u(1:j),u(2:j+1),'k.',r(1:j),r(2:j+1),'r.',lw,1);
  plot(vo,v,'k+',wo,w,'ro',lw,2); pause

end;


figure;
plot(t,u,'ro-',t,r,'bo-'); axis square
xlabel('Iteration Number, j','fontsize',16);
ylabel('Solution: v_j and w_j','fontsize',16);
title('Evolution of Frisch Map, different ICs','fontsize',16);
set (gcf,'color',[1 1 1]) % set the background color to white
set(gcf,'PaperUnits','normalized');set(gcf,'PaperPosition',[0 0 1 1])
print -dpng 'frisch_soln.png'
%!open frisch_soln.png

model = d0*2.^t; model=min(2,model);
semilogy(t,abs(r-u),'g.-',t,model,'linewidth',1.2); axis square
xlabel('Iteration Number, j','fontsize',16);
ylabel('Difference: |v_j-w_j|','fontsize',16);
title('Separation d_j=|v_j-w_j| and  2^jd_0 vs. j','fontsize',16);
set (gcf,'color',[1 1 1]) % set the background color to white
set(gcf,'PaperUnits','normalized');set(gcf,'PaperPosition',[0 0 1 1])
print -dpng 'frisch_growth.png'
%!open frisch_growth.png
