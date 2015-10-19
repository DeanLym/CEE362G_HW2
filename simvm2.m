function y = simvm2(N,Q,X)
%
%   y = simvm2(N,Q,X)
%
%   Generate N vectors of size ndata by 1 as realizations of 
%   a Gaussian vector with zero mean and covariance matrix 
%   Q(ndata, ndata)
%   If the third argument, matrix X(ndata,p), is given,
%   Q is interpreted as a generalized covariance matrix 
%   with X as the drift matrix.
%   This program uses eigenvalues decomposition and works with
%   (conditionally) semidefinite covariance matrices
%   
%   See also simv, simvm, simv2.

%   pkk 11/15/97
%
if nargin == 2
   [U,D] = eig(Q); d = diag(D); d = d.*(d>=0); 
   V = U*diag(sqrt(d));
   y = V*randn(max(size(Q)),N);
elseif nargin == 3
   [ndata,p] = size(X); [q,r] = qr(X); 
   q1 = q(1:ndata,1:p); q2 = q(1:ndata,p+1:ndata);
   Qa = q2'*Q*q2; Qa = (Qa+Qa')/2;
   [U,D] = eig(Qa); d = diag(D); d = d.*(d>=0); 
   V = U*diag(sqrt(d));
   ya = V*randn(max(ndata-p),N); y = q2*ya;
end