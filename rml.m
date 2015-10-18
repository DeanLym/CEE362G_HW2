function [theta, V, L] = rml(CT, theta, z, crit)
%[theta,V, L] = rml(CT, theta, z, crit)
%  Restricted Maximum Likelihood
%     for data covariance linear in np parameters theta 
%input:
%CT is a matrix, nz by nz by np, where C(:,:,1) is the derivative
%  with respect to parameter 1.
%theta is initial guess of the np parameters
%z is nz by 1 vector of data (generalized increments)
%crit is termination criterion (e.g., 0.001 or 0.0001)
%output:
%theta is vector with RML estimate
%V is its covariance matrix
%L is the value of the negative loglikelihood.

np = size(theta,1); nz = size(CT,1); PC = zeros(size(CT));
if size(CT,1)~=size(CT,2)
    error('CT must be nz by nz by np')
end
g = zeros(np,1); F = zeros(np);

k = 1; vcrit = 1/eps; L = 1/eps;
while (vcrit>crit)&(k<30)
    C = zeros(nz);
    for i=1:np
        C = C + CT(:,:,i)*theta(i);
    end
   iC = inv(C);
   for i = 1:np
      PC(:,:,i) = iC*CT(:,:,i);
   end
   for i = 1:np
      g(i) = 0.5*trace(PC(:,:,i))-0.5*z'*PC(:,:,i)*iC*z;
      for j = 1:np
         F(i,j) = 0.5*trace(PC(:,:,i)*PC(:,:,j));
      end
   end
   V = inv(F); 
   theta_new = theta - V*g;
   
   eigeC = eig(C); a = sum(log(eigeC)); L_up = 0.5*(a+z'*iC*z);
   if min(eigeC)<0
       warning('Loss of positive definiteness.')
   end
   vcrit = abs(L-L_up);
   theta = theta_new, L = L_up, g, k = k+1
   Q2 = z'*iC*z/nz
end

if k>=30
   warning('Algorithm may have not converged.')
end
Q2 = z'*iC*z/nz;
disp(['Test: Q2 is ', num2str(Q2), ' (should be 1).'])
