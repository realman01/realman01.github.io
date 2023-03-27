function R = moments_recursive(N,l,mu)
% 计算Legendre多项式L_0,L_1,...,L_N的矩
% lambda = (0 : 1 : N)';
R = zeros(mu+1,N+1);
for k = 0 : mu
    R(k+1,1) = gamma(k+1)/((1+l)^(k+1));
end
for n = 1 : N % 计算L_n的矩
    M = generateMatrix(l,mu,n);
    R(:,n+1) = M*(R(:,n));
end



function M = generateMatrix(l,mu,n)
% n >= 2
M = zeros(mu+1,mu+1);
for k = 0 : mu
    for j = 0 : k-1 % 对应log^{j}
        M(k+1,j+1) = (-1) * (gamma(k+1) / gamma(j+1)) * (2*n) / ((n+l+1)^(k-j+1));
%     M(k+1,j+1) = (-1) * compute_kj(k,j,n,l);
    end
    M(k+1,k+1) = (l-n+1)/(n+l+1);
end

function y = compute_kj(k,j,n,l)
y = 2*n/(n+l+1);
for i = j+1 : 1 : k
    y = y * (i/(n+l+1));
end
