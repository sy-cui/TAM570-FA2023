function Rp = r_periodic(N); % Periodic restriction matrix

N1=N+1;
Ih=speye(N1);
Rp=Ih(2:end,:); 
Rp(end,1)=1;
