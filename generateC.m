function c = generateC(lambda,beta)
% 关于 x^beta 权函数的积分
N = length(lambda);
c = zeros(N,1);
c(1) = 1/(lambda(1)+1+beta);
for n = 2 : N
    c(n) = c(n-1)*(lambda(n-1)-beta)/(-1-lambda(n)-beta);
end