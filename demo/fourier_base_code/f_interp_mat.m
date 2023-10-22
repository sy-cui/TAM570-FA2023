function [J] =  f_interp_mat(xo,xi);
%
%     Compute the Fourier interpolation matrix from xi to xo
%
%     ASSUMES xi points are UNIFORMLY SPACED, [x0,....,x_N],
%

      no = length(xo);
      ni = length(xi);
      N  = ni-1;        %% Number of Fourier basis fcts

      x0=xi(1); xN=xi(ni); L =xN-x0; dx=L/N;

      x = 2*pi*(xi(2:end)-x0)/L; %% On [0,2pi], start with x1
      z = 2*pi*(xo-x0)/L       ; %% On [0,2pi]

      Jh=ones ( N,N);
      Jf=ones (no,N);

      k=1;
      for kk=1:floor(N/2); 

          k=k+1; 


          Jh(:,k)=cos(kk*x);
          Jf(:,k)=cos(kk*z);

          k=k+1;
          if k<N+1;
            Jh(:,k)=sin(kk*x);
            Jf(:,k)=sin(kk*z);
          end;
      end;

J=Jf/Jh;
Avg=eye(N+1); Avg=Avg(2:end,:);
Avg(end,1)=0.5; Avg(end,end)=0.5;
J=J*Avg;

