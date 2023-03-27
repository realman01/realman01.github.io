function [xi,chi] = singular_gauss_quadrature(lamb,q,mu,N)
% guass-legendre quadrature for  P^{lamb,q}(log)^mu 
% N: number of quadrature nodes
l = (q+1)/lamb - 1;
N = N-1;
R = moments_recursive(2*N+1,l,mu);
abm = zeros(2*N+1,2);
for i = 0 : 2*N
    abm(i+1,1) = 1/2;
end
% abm(1,1) = 1;
for k = 1 : 2*N
    abm(k+1,2) = 1/4/(4-(k)^(-2));
end

mom = R(end,:)';
[alpha,beta]=my_chebyshev(N,mom,abm(:,1),abm(:,2));

A = zeros(length(alpha));
for i = 1 : length(alpha)
    A(i,i) = alpha(i);
end
for i = 1 : length(alpha)-1
    A(i,i+1) = sqrt(beta(i+1));
end
for i = 2 : length(alpha)
    A(i,i-1) = sqrt(beta(i));
end
[V,D] = eig(A);
t = diag(D); w = gamma(mu+1)/((l+1)^(mu+1)) * (V(1,:)' ./ sqrt(diag(V'*V))) .^2  ;


xi = t.^(1/lamb);
chi = xi.^(-q)  .* (-log(xi)).^(-mu) .* w .* (1/lamb)^(mu+1);

%% Test
% for k = 0 : 2* N+1
%     sum(xi.^(lamb*k+q) .* (-log(xi)).^(mu) .* chi) - gamma(mu+1)/((1+lamb*k+q).^(mu+1))
% end