function [x1, A1] = gaussQuad_weighted(lambda,beta,x,A,d,c)
N = length(lambda);
tol = 5*eps;
err = 1;
x1 = x;
A1 = A;
stepsize = 1;
iter = 0;
while err > tol
% % %     if mod(iter,5)==0
% % %     write('./log.txt', err);
% % %     end
% %     err
% %     x1 = reduce(x1);
% %     A1 = reduce(A1);
% %     [U,V,Y,Z] = generatesUVYZ(lambda,x1,A1);
% %     
% %     if iter > 1
% %        stepsize = stepsize * 1;
% % %         beta = 0.8;
% %     end
% %     S = (V*(U\Y)-Z);
% %     s = (d-V*(U\c));
% %     a = S\s;
% %     direction = beta*direction + (1-beta)*a;
% %     x2 = x1 - stepsize * direction;
% %     A2 =  U\(Y*(x1-x2)+c);
%% v1
% err
[U,V,Y,Z] = generatesUVYZ(lambda+beta/2,x1,A1);
B = zeros(N,N);
B(1:N/2,1:N/2) = U;
B(N/2+1:N,1:N/2) = V;
B(N/2+1:N,N/2+1:N) = (-beta/2) .* V + Z;
B(1:N/2,N/2+1:N) = (-beta/2) .* U + Y;
rk = ([U;V]*diag(1./(x1.^(beta/2)))*A1 - [c;d]); %%%%%%%%% 修改beta/2 --> 1+beta/2
p = B \(-rk);

if iter > 10
        stepsize = stepsize * 0.1;
end
% temp = [A1;x1] + stepsize*p;
% A2 = temp(1:N/2,1);
% x2 = temp(N/2+1:end);

delta_A = p(1:N/2,1);
delta_x = p(N/2+1:N,end);
A2 = A1 + stepsize*delta_A.*x1.^(beta/2); %%%%%%%%% 修改beta/2 --> 1+beta/2
x2 = x1 + stepsize * ((x1.^(1+beta/2))./A1) .* delta_x; %%%%% 修改1+beta/2 --> 2+beta/2

    
%%
    
    err = norm(x1-x2);
%     err = max(abs(rk));
%     [U,V,~,~] = generatesUVYZ(lambda,x1,A1);
%     err = max(abs([U;V]*A2 - [c;d]));
    x1 = x2;
    A1 = A2;
    iter = iter + 1;
end