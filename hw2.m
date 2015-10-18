

function [] = hw2()
    clear all;
    close all;
    data = load('2015Assign1_1.txt');
    t = data(:,1);
    y = data(:,2);
    % number of unknowns and number of observations
    m = 400; n = length(y);
    % get the integral matrix
    H = shaw(m); 
    X = ones(m,1);
    % Q0 and R0
    Q0 = eye(m); R0 = eye(n);
    % PHI
    PHI = H*X;
    % Transformation
    T = null(PHI')';
    z = T*y; nz = length(z);
    CT = zeros(nz,nz,2);
    CT(:,:,1) = T*H*Q0*H'*T';
    CT(:,:,2) = T*T';
    [theta_rml,Vp, L] = rml(CT, [5;5], z, 0.0001)
    % the negative loglikelihood function L(theta)
    function res = GetL(theta1,theta2)
        PSI = H*theta1*Q0*H'+theta2*R0;
        TPSI = T*PSI*T';
        log_det_TPSI = sum(log(eig(TPSI)));
        res = 0.5*log_det_TPSI+0.5*z'*inv(TPSI)*z;
    end
    
%% Question 1
n_theta = 100;
theta1 = linspace(18,20,n_theta);
theta2 = linspace(0.5e-7,2e-4,n_theta);
L = zeros(n_theta);
X = zeros(n_theta,1);
Y = zeros(n_theta,1);
for i=1:n_theta
   for j = 1:n_theta
      L(i,j) = GetL(theta1(i),theta2(j));
      X(i) = theta1(i);
      Y(j) = theta2(j);
   end
end
stop = 1;

figure;
contour(L);
figure;
surf(X,Y,L,L);
end



function [] = q1()
    




end


