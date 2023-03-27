function [x1,A1] = continuation_weighted(lambda,stepsize,beta,x1,A1,polys)
% lambda(1) must be 0 !!!!
% Consider continuation applying on Lambda
N = length(lambda);
%%
if nargin == 3
%     [x1,A1] = Legendre_Nodes(N/2,0,1);
    [t,w] = Jacobi_Nodes(0,beta,N/2);
    x1 = (t+1)/2;
    A1 = w/(2^(1+beta));
    polys = (0:N-1)';
    clear t w
end
%%
% if nargin == 2
%     lamb = lambda(3);
%     [x1,A1] = Jacobi_Muntz_Nodes(lamb,0,N/2);
%     polys = (0:N-1)'*lamb;
% end
%%
% if nargin == 2
%     [x1,A1] = temp1();
%     polys = zeros(30,1);
%     for i = 1 : 60
%         polys(i) = floor((i-1)/2);
%     end
% end
%%
% spline interpolation
% x = [0 1];
% y = [1 0];
% cs = spline(x,[0 y 0]);
% alpha = 0 : stepsize : 1;
% if alpha(end) ~= 1
%     alpha(end+1) = 1;
% end
% alpha = ppval(cs,alpha);

%%
alpha = 0 : stepsize : 1;
if alpha(end) ~= 1
    alpha(end+1) = 1;
end
alpha = alpha(end:-1:1);
%%

% h=waitbar(0,'please wait');

for i = 1 : length(alpha)
    t = generateC((alpha(i)).*polys+(1-alpha(i)).*lambda + beta/2,beta/2);
    d = t(N/2+1:N,1);
    c = t(1:N/2,1);
    [x1, A1] = gaussQuad_weighted((alpha(i)).*polys+(1-alpha(i)).*lambda,beta,x1,A1,d,c);
%     if flag == 1
%
%     elseif flag == 0
%         if alpha(i) < alpha(i-1)
%             alpha(i) = alpha(i) + stepsize^2;
%         else
%             disp("ERROR");
%             return
%         end
%     end
    
    str=['Continuation for Gauss Quadrature...',num2str(i/length(alpha)*100),'%'];
%     waitbar(i/length(alpha),h,str)
    disp(str)
end

% delete(h);