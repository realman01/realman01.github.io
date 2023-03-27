function [val, diff] = muntz_legendre3(lambda,x)
p.x = x; 
p.w = -log(x);
p.lambda = lambda;

% p.sigma = min(lambda) + min(1./log(x),-1/5);

p.sigma = zeros(size(x));
for i = 1 : length(x)
    [s,~,~,~] = search(lambda,-log(x(i)));
    p.sigma(i) = min(lambda) + max(s,1)/log(x(i));

end

%%
% p.m >= number of Gauss quadrature nodes
p.m = 80;%floor( sqrt( 4*(2*max(lambda)+1)*max(p.w)+1 )/2 + 1/2) +1;
p.m = ceil(abs(p.m))+1;
p.a = p.m*pi;
p.glen = 80;%max(20,ceil(length(lambda))); % # gauss legendre nodes
p.glan = 80;%max(20,ceil(length(lambda))); % # gauss laguerre nodes
[p.gle_t, p.gle_w] = Legendre_Nodes(p.glen,0,1); 
[p.gla_t, p.gla_w] = Laguerre_Nodes(p.glan, 0);

le = legendre_part(p);
la = laguerre_part(p);
val = diag(x.^(p.sigma)) * imag(le+la)./pi; %%%%%%%%%%%% 修改p.sigma+1

if nargout > 1
%    disp('Caculating derivatives...');
   diff = zeros(length(x),length(lambda));
   diff(:,1) = lambda(1) .* x.^(lambda(1));
   for k = 2 : length(lambda)
      diff(:,k) = diff(:,k-1) + ( lambda(k) * val(:,k) + (1 + lambda(k-1))* val(:,k-1));
   end
%    diff = diag(x) \ diff;
end

%% Compute Legendre Part
function le = legendre_part(p)
le = zeros(length(p.x),length(p.lambda));
H = zeros(length(p.x),p.glen,p.m);
% val(:,1) = 
for i = 1 : length(p.lambda)
    H = generate_H(p,i,H);
%     for k = 1 : p.m
%         le(:,i) = le(:,i) + H(:,:,k)*p.gle_w;
%     end
    le(:,i) = sum(H,3)*p.gle_w;
end
le = pi*le;



function H2 = generate_H(p,n,H1)
% generate 3d matrix H
% n>2
if n>length(p.lambda),return; end
t = p.gle_t;
w = p.w;
w = reshape(w,length(w),1); % x column vector
t = reshape(t,1,length(t)); % t row vector
if n==1
    H2 = zeros(length(p.w),p.glen,p.m);
    for k = 1 : p.m
        H2(:,:,k) = 1./(pi*(t+k-1) + 1i*(p.sigma-p.lambda(1)).*w) *...
            diag(exp(1i*pi*(t+k-1)'));
    end
    return;
end

W = zeros(length(p.w),p.glen,p.m);
for k = 1 : p.m
    W(:,:,k) = (pi*(t+k-1) + 1i*(p.sigma + p.lambda(n-1) + 1).*w) ./ (pi*(t+k-1) + 1i*(p.sigma-p.lambda(n)).*w);
end

H2 = H1.*W;


%% Compute Laguerre part
function la = laguerre_part(p)
la = zeros(length(p.x),length(p.lambda));
G = zeros(length(p.w),p.glan);
for i = 1 : length(p.lambda)
    G = generate_G(p,i,G);
    la(:,i) = G * (exp(-p.gla_t) .* p.gla_w);
end
la = 1i*(-1)^p.m*la;

function G2 = generate_G(p,n,G1)
% generate 3d matrix H
% n>2

if n>length(p.lambda),return; end
t = p.gla_t;
w = p.w;
w = reshape(w,length(w),1); % x column vector
t = reshape(t,1,length(t)); % t row vector
if n==1
    G2 = 1./(p.a+1i*t + 1i*(p.sigma-p.lambda(1)).*w);
    return;
end

W = (p.a+1i*t + 1i*(p.sigma + p.lambda(n-1) + 1).*w) ./ (p.a+1i*t + 1i*(p.sigma-p.lambda(n)).*w);

G2 = G1.*W;