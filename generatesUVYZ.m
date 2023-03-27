function [U,V,Y,Z] = generatesUVYZ(lambda,x,A)
N = length(lambda);
if mod(N,2) ~= 0
    disp('Lambda Length Not Even!!');
end
[val, diff] = muntz_legendre3(lambda,x);
U = val(:,1:N/2)';
V = val(:,N/2+1:N)';
Y = diff(:,1:N/2)';% * diag(A);
Z = diff(:,N/2+1:N)';% * diag(A);