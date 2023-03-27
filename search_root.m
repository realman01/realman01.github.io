function [s,fval,exitflag,output] = search_root(lambda,w)
y = @(theta) myfun(lambda,w,theta);
% [s,fval,exitflag,output] = fminbnd(y,0,20);
% [s,fval,exitflag,output] = fminunc(y,pi);
[s,fval,exitflag,output] = fminsearch(y,pi);

function y = myfun(lambda,w,theta)
% y = 0;
% for i = 1 : length(lambda) - 1
%     y = y + log(abs(sigma+lambda(i)+1)) - log(abs(sigma-lambda(i)));
% end
% 
% y = y - log(abs(sigma-lambda(i)));
% y = y - log(w);
% y = y + sigma.*(-w);

%%
y = ones(size(w));
for i = 1 : length(lambda)-1
    y = y .* (theta - (lambda(i)+min(lambda)+1)*w) ./ (theta+(min(lambda)+lambda(i))*w);
end    
y = y ./ abs(theta + (min(lambda)+lambda(end)).*w);
y = abs(y) .* exp(sqrt(w)) +  exp(-min(lambda)*w) .* exp(theta) ./ abs(sqrt(theta));