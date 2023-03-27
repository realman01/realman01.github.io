% n个节点的Gauss积分公式，区间(-1,1)
function [t,w] = Jacobi_Nodes(alpha,beta,n)
if alpha == -1/2 && beta == -1/2 
    [t,w] = Chebyshev_Nodes(n,-1,1);
    return; 
end

for j = 1 : 1 : n-1
    a(j) = (beta^2-alpha^2) ./ ((2*j + alpha + beta) * (2*j + alpha + beta + 2));
    b(j) = 4*j*(j+alpha)*(j+beta)*(j+alpha+beta) ./ ((2*j+alpha+beta-1) * (2*j+alpha+beta+1) * (2*j+alpha+beta)^2);
end
temp = beta-alpha;
temp = temp /(alpha+beta+2);
a = [temp,a];
b = b'; a = a';
A = diag(a) + diag(sqrt(b),1) + diag(sqrt(b),-1);
[V,D] = eig(A);
t = diag(D); w = V(1,:).^2;
w = w .* (2^(alpha + beta + 1) * gamma(alpha+1) * gamma(beta + 1)) ./ gamma(alpha+beta+2);
w = w';