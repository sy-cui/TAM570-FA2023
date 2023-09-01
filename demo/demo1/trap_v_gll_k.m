hdr

% \int_0^1 sin(k pi x) dx = (1-cos(k pi ))/(k pi)

I_exact = (1-cos(k*pi))/(k*pi);

hold off
n=200+8*ceil(k); x=(0:n)/n; f=sin(k*pi*x); plot(x,0*x,'k-',x,f,'k-',lw,1.1);

%  Trapezoidal rule test

for n=2:200;

   x=(0:n)/n; f=sin(k*pi*x);
   I_trap(n-1) = ( sum(f(2:n)) + .5*(f(1)+f(n+1)) ) / n;
   e_trap(n-1) = abs(I_exact - I_trap(n-1));
   n_trap(n-1) = n;

   if n==20;
      hold on
      plot(x,f,'r+',lw,2);
      title('Trapezoidal Rule, N=20',fs,16);
      xlabel('x',fs,16);
      ylabel('f(x)',fs,16);
   end;

end;


%  GLL Test

for N=2:100;
   [z,w] = zwgll(N); x = .5*(z+1);
   f = sin(k*pi*x);
   I_gll(N-1)=.5*sum(w.*f);
   e_gll(N-1)=abs(I_exact-I_gll(N-1))+eps;
   n_gll(N-1)=N;

   if N==20;
      pause(.1);
      hold on
      plot(x,f,'bx',lw,2);
      title('Trapezooidal and GLL Quadrature, N=20',fs,16);
      xlabel('x',fs,16); ylabel('f(x)',fs,16);
      hold off
   end;

end;
pause;

figure
loglog(n_trap,e_trap,'ro-',lw,2); axis square;
y=1.e-14; x=200; st=['k = ' num2str(k)]; text(x,y,st,fs,14)
title('Trapezoidal Rule Convergence: \int sin(k \pi x)',fs,16)
xlabel('N',fs,16); ylabel('Integration Error',fs,16); 
legend('Composite Trapezoidal Rule','location','southwest')
pause

figure
loglog(n_gll,e_gll,'bo-',n_trap,e_trap,'ro-',lw,2); axis square;
y=1.e-12; x=2; st=['k = ' num2str(k)]; text(x,y,st,fs,14)
title('GLL and Trapezoidal Rule Convergence: \int sin(k \pi x)',fs,16)
xlabel('N',fs,16); ylabel('Integration Error',fs,16); 
legend('Gauss-Lobatto-Legendre','Composite Trapezoidal Rule','location','southwest')
pause

figure
semilogy(n_gll,e_gll,'bo-',n_trap,e_trap,'ro-',lw,2); axis square;
title('GLL and Trapezoidal Rule Convergence: \int sin(k \pi x)',fs,16)
xlabel('N',fs,16); ylabel('Integration Error',fs,16); 
legend('Gauss-Lobatto-Legendre','Composite Trapezoidal Rule','location','northeast')
y=1.e-12; x=100; st=['k = ' num2str(k)]; text(x,y,st,fs,14)

