%-------------------------------------------------------
lw  =   'linewidth';
fs  =   'fontsize';
intp =  'interpreter';
ltx  =  'latex';
format compact;
%%% Plotting defs

% \int_0^1 sin(k pi x) dx = (1-cos(k pi ))/(k pi)
k = 3.1
I_exact = (1-cos(k*pi))/(k*pi);
hold off

n=200+8*ceil(k); x=(0:n)/n; f=sin(k*pi*x); plot(x,f,'k-',lw,1.1);

for n=2:500;                           %  TRAPEZOIDAL RULE TEST
   x=(0:n)/n; f=sin(k*pi*x);
   I_trap(n-1) = ( sum(f(2:n)) + .5*(f(1)+f(n+1)) ) / n;
   e_trap(n-1) = abs(I_exact - I_trap(n-1));
   n_trap(n-1) = n;
   if n==20; hold on;
      plot(x,f,'r+',lw,2);
      title('Trapezoidal Rule, N=20','FontSize',16);
      xlabel('x','FontSize',16);
      ylabel('f(x)','FontSize',16);
      pause;
    end; 
end;

for N=2:400;
   [z,w] = zwgll(N); x = .5*(z+1);
   f = sin(k*pi*x);
   I_gll(N-1)=.5*sum(w.*f);
   e_gll(N-1)=abs(I_exact-I_gll(N-1))+eps;
   n_gll(N-1)=N;
   if N==20; hold on;
      plot(x,f,'bx',lw,2);
      title('GLL Quadrature, N=20','FontSize',16);
      xlabel('x','FontSize',16);
      ylabel('f(x)','FontSize',16);
      hold off; pause;
end; end;
loglog(n_trap,e_trap,'ro-'); axis square;
title('Trapezoidal Rule Convergence: \int sin(k \pi x)','FontSize',16)
xlabel('N','FontSize',16); ylabel('Integration Error','FontSize',16);
legend('Composite Trapezoidal Rule','location','southwest'); pause
loglog(n_gll,e_gll,'bo-',n_trap,e_trap,'ro-'); axis square;
title('GLL and Trapezoidal Rule Convergence: \int sin(k \pi x)','FontSize',16)
xlabel('N','FontSize',16); ylabel('Integration Error','FontSize',16);
legend('Gauss-Lobatto-Legendre','Composite Trapezoidal Rule','location','southwest'); pause
semilogy(n_gll,e_gll,'bo-',n_trap,e_trap,'ro-'); axis square;
title('GLL and Trapezoidal Rule Convergence: \int sin(k \pi x)','FontSize',16)
xlabel('N','FontSize',16); ylabel('Integration Error','FontSize',16);
legend('Gauss-Lobatto-Legendre','Composite Trapezoidal Rule','location','southwest'); pause

%-------------------------------------------------------