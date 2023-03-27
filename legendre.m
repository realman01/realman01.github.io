function val = legendre(n,lambda,x)
% val = legendre(n,lambda,x) compute nth order muntz-legendre poly's value at x
% Input: n is a integer s.t. 1 <= n <= length(lambda) 
%        lambda is a column vector
%        x is a column vector
% Output: 
%        val: a vector as same size as x
%% Empirical parameters
p.m = 60;
p.h = pi;
p.a=p.m*p.h;
p.num_nodes = 30;
%%
p.n=n;
p.sigma = zeros(size(x));
for i = 1 : length(x)
    [s,fval,exitflag,output] = search_root(lambda(1:n),-log(x(i)));
    p.sigma(i) = min(lambda) + max(s,1)/log(x(i));
end
l1 = L1(p,lambda,x);
l2 = L2(p,lambda,x);
val = x.^(p.sigma) .* imag(l1+l2)./pi;

function s= L1(p,lambda,x)
s = zeros(size(x));
[nodes,weights] = Legendre_Nodes(p.num_nodes,0,1);
for k=1:p.m
    quad = 0;
    for j = 1:length(nodes)
        quad = quad + f(p,lambda,p.h.*(nodes(j)+k-1),x).*exp(1i*p.h*(nodes(j)+k-1)).*weights(j);
    end
    s = s + quad;
end
s = s*p.h;
    
function s = L2(p,lambda,x)
[nodes, weights] = Laguerre_Nodes(p.num_nodes, 0);
s = 0;
for j = 1 : length(nodes)
    s = s + phi(p,lambda,nodes(j),x).*exp(-nodes(j)).*weights(j);
end
s = exp(1i*p.a) * s;

function s = phi(p,lambda,y,x)
s = 1i*f(p,lambda,p.a+1i*y,x);

function s = f(p,lambda,t,x)
% t=reshape(t,length(t),1);
w = -log(x);
s = ones(size(t));
for j = 1:p.n-1
    s = s.*(t+1i.*(p.sigma+lambda(j)'+1).*w) ./ (t+1i*(p.sigma-lambda(j)).*w);
end
s = s ./ (t+1i*(p.sigma-lambda(p.n)).*w);

