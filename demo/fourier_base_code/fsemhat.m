      function[Ah,Bh,Ch,Dh,z,w] = fsemhat(N)
%
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a trigonometric polynomial
%     of degree N
%

      [z,w] = zwuni(N);

      Bh    = diag(w);
      Dh    = f_deriv_mat(z);

      Ah    = Dh'*Bh*Dh;
      Ch    = Bh*Dh;
