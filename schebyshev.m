function [ab,normsq]=schebyshev(dig,N,mom,abm)
%SCHEBYSHEV Modified Chebyshev algorithm.
%   [AB,NORMSQ]=SCHEBYSHEV(DIG,N,MOM,ABM) uses the modified
%   Chebyshev algorithm to generate the Nx2 array AB of the first 
%   N recurrence coefficients for polynomials orthogonal with
%   respect to a weight function specified by the 1x(2N) array
%   MOM of its first 2N modified moments relative to monic
%   polynomials defined by the (2N-1,2) array ABM of their
%   recurrence coefficients. The squared norms of these orthogonal
%   polynomials are placed into the 1xN array NORMSQ. The call
%   [AB,NORMSQ]=SCHEBYSHEV(DIG,N,MOM) does the same, but using the
%   classical Chebyshev algorithm. If N is larger than the sizes
%   of MOM and ABM warrant, then N is reduced accordingly.

syms sig
digits(dig); 
if N<=0, error('N out of range'), end
if N>size(mom,2)/2, N=size(mom,2)/2; end
if nargin<4, abm=zeros(2*N-1,2); end
if N>(size(abm,1)+1)/2; N=(size(abm,1)+1)/2; end
smom=vpa(mom,dig); sabm=vpa(abm,dig);
ab(1,1)=sabm(1,1)+smom(2)/smom(1); ab(1,2)=smom(1);
if N==1, normsq(1)=smom(1); return, end
sig(1,1:2*N)=0; sig(2,:)=smom(1:2*N);
for n=3:N+1
  for m=n-1:2*N-n+2
    sig(n,m)=sig(n-1,m+1)-(ab(n-2,1)-sabm(m,1))*sig(n-1,m) ...
      -ab(n-2,2)*sig(n-2,m)+sabm(m,2)*sig(n-1,m-1);
  end
  ab(n-1,1)=sabm(n-1,1)+sig(n,n)/sig(n,n-1)-sig(n-1,n-1)/ ...
    sig(n-1,n-2);
  ab(n-1,2)=sig(n,n-1)/sig(n-1,n-2);
end
for n=1:N
  normsq(n)=sig(n+1,n); 
end

%% 
% Example
% find the nodes and weights in Gauss-Laguerre quadrature

% A =diag(ab(:,1)) +  diag(sqrt(ab(2:end,2)),1) + diag(sqrt(ab(2:end,2)),-1);
% [V,D] = eig(A);
% t = diag(D); w = V(1,:).^2;
% w = w .* 2;
% w = w';

