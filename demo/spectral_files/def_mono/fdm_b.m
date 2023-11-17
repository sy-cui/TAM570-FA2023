function [Zb]=fdm_b(Fb,Rr,Rs,Sx,Sy,Di);

Zb=Rr'*(Sx*(Di.*(Sx'*(Rr*Fb*Rs')*Sy))*Sy')*Rs; 

