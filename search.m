function [s,fval,exitflag,output] = search(lambda,w)
%% search(lambda,w) find the optimal sigma in computing Müntz Legendre polynomials.
% Input: lambda is a column vector
%        w is a real number
% Output: s is the optimal sigma. fval, exitflag and output are the parameters corresponding to fminsearch.m        
y = @(theta) myfun(lambda,w,theta);
[s,fval,exitflag,output] = fminsearch(y,pi);

function y = myfun(lambda,w,theta)
y = ones(size(w));
for i = 1 : length(lambda)-1
    y = y .* (theta - (lambda(i)+min(lambda)+1)*w) ./ (theta+(min(lambda)+lambda(i))*w);
end    
y = y ./ abs(theta + (min(lambda)+lambda(end)).*w);
y = abs(y) .* exp(sqrt(w)) +  exp(-min(lambda)*w) .* exp(theta) ./ abs(sqrt(theta));
