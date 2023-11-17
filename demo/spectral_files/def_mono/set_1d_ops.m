function [S,D]=set_1d_ops(Ah,Bh,R,L);

Li=2/L; L2=L/2;

A=Li*R*Ah*R';
B=L2*R*Bh*R'; 

[S,D]=gen_eig_ortho(A,B);


% S=R'*S;   %  Prolongate to save on subsequent multiplies


