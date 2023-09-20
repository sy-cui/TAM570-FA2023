function [S,Lam]=gen_eig_decomp(A,B);

%
%  Return matrix of eigenvectors,   S   = [s1 ... sn],
%     and matrix of eigenvalues,    Lam = diag(lam_k).
%
%  The columns of S are the eigenvectors, sk, k=1,...,n
%  Lam is defined as a sparse matrix.
%
%
%  The eigenvalues are sorted small to large
%  The eigenvectors are orthonormalized w.r.t. the B inner-product:
%
%  sj^T B sk = delta_jk
%


%   Below is new code provided by James Lottes
A=full(A+A')/2;
B=full(B+B')/2;
% [S,Lam] = eig(A,B,'chol');  %% This doesn't seem to work for Octave
[S,Lam] = eig(A,B);
lambda = diag(Lam);
[lambda, i] = sort(lambda);
S = S(:, i);
S = S / chol(S' * B * S);
S = S / chol(S' * B * S);  % 2nd orthogonalization may not be necessary

Lam=diag(lambda); Lam=sparse(Lam);


