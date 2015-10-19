function y = simv(Q,X)
%
%   y = simv(Q,X)
%
%   Generate a vector of size ndata by 1 as a realization of 
%   a Gaussian vector with zero mean and covariance matrix 
%   Q(ndata, ndata)
%   If the second argument, matrix X(ndata,p), is given,
%   Q is interpreted as a generalized covariance matrix 
%   with X as the drift matrix.
%   This program uses Cholesky decomposition and requires
%   (conditionally) positive definite covariance matrices.
%
%   See also simvm, simv2, simvm2.

%   pkk 11/15/97
%
if nargin == 1
   V=(chol(Q))';
   y = V*randn(max(size(Q)),1);
elseif nargin ==2
   [ndata,p] = size(X); [q,r] = qr(X); 
   q1 = q(1:ndata,1:p); q2 = q(1:ndata,p+1:ndata);
   Qa = q2'*Q*q2; V=(chol(Qa))';
   ya = V*randn(max(ndata-p),1); y = q2*ya;
end