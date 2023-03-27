function [val, diff] = muntz_legendre(lambda,x)
% muntz_legendre(lambda,x) evaluate Müntz Legendre polynomials and their first-order derivatives for all lambda.
% Input:
%   x: N-vector, the nodes in which need to be estimated
%   lambda: n-vector, the muntz poly orders
%   p: parameters structure
%       p.num_nodes
%       p.m
% Return:
%   val: N by n matrix, ith column cotains m values of ith Muntz Legendre poly at x
%%
val = zeros(length(x),length(lambda));
% val(:,1) = ones(length(x),1);%
for k = 1 : length(lambda)
   val(:,k) = legendre(k,lambda,x);
end

if nargout > 1
%   disp('Caculating derivatives...');
   diff = zeros(length(x),length(lambda));
   diff(:,1) = lambda(1) .* x.^(lambda(1)-1);
   for k = 2 : length(lambda)
      diff(:,k) = diff(:,k-1) + ( lambda(k) * val(:,k) + (1 + lambda(k-1))* val(:,k-1))./x;
   end
end
