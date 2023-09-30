addpath('/Users/sycui/Desktop/course_work/Fall_2023/TAM570/lib')

premeable;
ns = 21:60;

errors = ns;
for i = 1:length(ns);
    N = ns(i);
    Lx=1;

    [Ah,Bh,Ch,Dh,z,w] =  semhat(N);
    Ih=eye(N+1);

    Bb = (Lx/2)*Bh;
    Ab = (2/Lx)*Ah;

    R=Ih(2:N,:);

    A=R*Ab*R';
    B=R*Bb*R';

    [V,D]=eig(A,B);

    lam=diag(D);
    lam=sort(lam);

    k = [1:N-1]';
    lt = pi*pi*(k.*k);   % Exact

    % plot(k,lt,'bo',lw,2,k,lam,'ro',lw,2);
    error_curr = abs(lt - lam) ./ lt; 
    errors(i) = error_curr(20);

    

end;

figure
semilogy(ns, errors, 'ro', lw, 2);
xlabel("N")
ylabel("Error")