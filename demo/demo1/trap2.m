clear all
hdr

% \int_0^1 sin(k pi x) dx = (1-cos(k pi ))/(k pi)



%  Trapezoidal rule test

for kk=1:3;
    if kk==1; clr='r-'; k= 3.1; end;
    if kk==2; clr='b-'; k=31.1; end;
    if kk==3; clr='k-'; k=61.1; end;

I_exact = (1-cos(k*pi))/(k*pi);

%hold off
%n=200+8*ceil(k); x=(0:n)/n; f=sin(k*pi*x); 
% plot(x,0*x,'k-',x,f,'k-',lw,1.1);

for N=2:200;

   x=(0:N)/N; f=sin(k*pi*x);
   I_trap(N-1) = ( sum(f(2:N)) + .5*(f(1)+f(N+1)) ) / N;
   e_trap(N-1) = abs(I_exact - I_trap(N-1));
   n_trap(N-1) = N;

%  if N==20;
%     hold on
%     plot(x,f,'r+',lw,2);
%     title('Trapezoidal Rule, N=20',fs,16);
%     xlabel('x',fs,16);
%     ylabel('f(x)',fs,16);
%  end;

end;


%  GLL Test

figure(2);
for N=2:200;
   [z,w] = zwgll(N); x = .5*(z+1);
   f = sin(k*pi*x);
   I_gll(N-1)=.5*sum(w.*f);
   e_gll(N-1)=abs(I_exact-I_gll(N-1))+eps;
   n_gll(N-1)=N;

%  if N==20;
%     pause(.1);
%     hold on
%     plot(x,f,'bx',lw,2);
%     title('GLL Quadrature, N=20',fs,16);
%     xlabel('x',fs,16); ylabel('f(x)',fs,16);
%     hold off
%  end;

end;

figure(1);
semilogy(n_trap,e_trap,clr,lw,2); axis square; axis([1 200 1e-5 1]);
y=1.e-14; x=200; st=['k = ' num2str(k)]; %text(x,y,st,fs,14)
title(['Trapezoidal Rule: k=3.1, 31.1, 61.1'],fs,16)
xlabel('N',fs,16); ylabel('Integration Error',fs,16); 
% legend('Composite Trapezoidal Rule','location','southwest')
hold on;
semilogy(n_trap,.001+0*n_trap,'k--',lw,1); 
pause(.03)


figure(2)
semilogy(n_gll,e_gll,clr,lw,2); axis square; axis([1 200 1e-16 1]);
y=1.e-12; x=2; st=['k = ' num2str(k)]; %text(x,y,st,fs,14)
%title(['GLL: ' num2str(k)],fs,16)
title('Gauss Lobatto Legendre: k=3.1, 31.1, 61.1',fs,16)
xlabel('N',fs,16); ylabel('Integration Error',fs,16); 
% legend('Gauss-Lobatto-Legendre','location','southwest')
hold on;
semilogy(n_gll,.001+0*n_gll,'k--',lw,1); 
pause(.03)

end;

