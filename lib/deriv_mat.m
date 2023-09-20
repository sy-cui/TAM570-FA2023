      function[d] =  deriv_mat(x)
%
%     Compute the interpolation matrix from x to x
%
      
      n = length(x);
      a = ones(n,1);
      for i=1:n;
        for j=1:(i-1);  a(i)=a(i)*(x(i)-x(j)); end;
        for j=(i+1):n;  a(i)=a(i)*(x(i)-x(j)); end;
      end;
      a=1./a; % These are the alpha_i's

      d = zeros(n,n);
      for j=1:n; for i=1:n;
         if i~=j; d(i,j) = a(j)/( a(i)*(x(i)-x(j)));end;
      end;end;

      for i=1:n; d(i,i)=0; d(i,i)=-sum(d(i,:)); end;
