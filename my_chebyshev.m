function [alpha,beta]=my_chebyshev(N,mom,a,b)
for l = 0 : 2*N+1
    A(1,l+1) = 0;
    A(2,l+1) = mom(l+1);
end
alpha = zeros(N+1,1); beta = zeros(N+1,1);
alpha(1) = a(1) + (1/2)*mom(2)/mom(1);
for k = 1 : N 
    for l = k : 2*N+1-k
        A(k+2,l+1) = (A(k+1,l+2)*(l+1)/(4*l+2) + A(k+1,l) * b(l+1) * (4*l-2)/l) +...
            ((a(l+1)-alpha(k))*A(k+1,l+1) - beta(k) * A(k,l+1));
    end
    alpha(k+1) = a(k+1) - (k/(4*k-2))*(A(k+1,k+1)/A(k+1,k)) + (A(k+2,k+2)/A(k+2,k+1)) * ((k+1)/(4*k+2));
    beta(k+1) = A(k+2,k+1)/A(k+1,k) * (k/(4*k-2));
end
beta(1) = A(2,1);
