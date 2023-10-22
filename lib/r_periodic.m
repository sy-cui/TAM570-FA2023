function Rp = r_periodic(N); % Periodic restriction matrix

Rp=speye(N+1)(2:end,:);
Rp(end,1)=1;
