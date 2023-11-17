function[D] =  f_deriv_mat(x)
hdr

%
%     x = [x0 ... xN]
%
%     ASSUMES x points are UNIFORMLY SPACED, [x0,....,x_N],
%

      x0 = x(1); xN=x(end);

      L  = xN-x0;

      n = length(x);
      N = n-1;        %% Number of Fourier basis fcts

      x0=x(1); xN=x(n); L =xN-x0; dx=L/N;

      x = 2*pi*(x-x0)/L; %% On [0,2pi], start with x1
      z = x(2:end);

      D =zeros(n,N);
      J =ones (N,N);

      k=1;
      for kk=1:floor(N/2); 
          k=k+1;           
            J(:,k) =     cos(kk*z);
            D(:,k) = -kk*sin(kk*x);
          k=k+1; if k<N+1; 
            J(:,k) =     sin(kk*z);
            D(:,k) =  kk*cos(kk*x); 
          end;
      end;

Avg=eye(n); Avg=Avg(2:end,:);
Avg(end,1)=0.5; Avg(end,end)=0.5;

D = (2*pi/L)*(D*inv(J));
D = D*Avg;

