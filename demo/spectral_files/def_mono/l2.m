function W = l2(Ub,Rr,Rs,Bb)

U=Rr*(Bb.*Ub)*Rs';

B=Rr*Bb*Rs'; Bi=1./B;

W=Rr'*(Bi.*U)*Rs;

