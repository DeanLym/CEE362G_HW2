function [s,V,LAMBDA,MU] = GenLinInv(y,H,R,X,Q)
%  function [s,V,LAMBDA,MU] = GenLinInv(y,H,R,X,Q)
%  
%       General Linear Estimation.
%       Calculates best estimates and covariance matrix of an
%       unknown vector using linear inverse problem
%       also computes LAMBDA and MU
%
%  y - vector of measurements
%  H - observation matrix
%  R - Covariance matrix for measurement residual
%  X - Drift matrix for s
%  Q - Covariance matrix for unknowns
%  s - Best estimates (posterior mean) for unknowns
%  V - Error (posterior) covariance matrix for unknowns. Diagonal of
%        this matrix contains the mean square estimation errors.

[m,p]=size(X);
[n,m1] = size(H);
[n1,n2] = size(R);
[m2,m3] = size(Q);

if var([m,m1,m2,m3]>0)
    error('check dimension m')
elseif var([n,n1,n2]>0)
    error('check dimension n')
end


rH = rank(H); T = [];
if rH<n
    [UH,SH,VH] = svd(H);
    T = UH(:,1:rH)';
    n = rH;
    H = T*H; 
    R = T*R*T';
    y = T*y;
end

PHI = H*X; QHT = Q*H'; HQHT = H*QHT;
PSI = HQHT+R;

AA = [PSI, PHI; PHI', zeros(p)]; AA=(AA+AA')/2;
bb = [QHT';X'];

SOL = AA\bb; 
LAMBDA = (SOL(1:n,:))'; 
MU = SOL(n+1:n+p,:);

s = LAMBDA*y; 
V = -X*MU+Q-QHT*LAMBDA';

if ~isempty(T)
    LAMBDA = LAMBDA*T;
end    






