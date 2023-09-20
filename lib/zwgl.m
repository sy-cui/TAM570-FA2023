% GAUSS Gauss quadrature rule.
%
%    Given a weight function w encoded by the nx2 array ab of the 
%    first n recurrence coefficients for the associated orthogonal
%    polynomials, the first column of ab containing the n alpha-
%    coefficients and the second column the n beta-coefficients, 
%    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
%    the n-point Gauss quadrature rule for the weight function w.
%    The nodes, in increasing order, are stored in the first 
%    column, the n corresponding weights in the second column, of
%    the nx2 array xw.
%
function [z,w]=zwgl(N)
J=eye(N);
for n=2:N
  J(n,n-1)=1;
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];
[z,w]=[xw(:,1) xw(:,2)];
