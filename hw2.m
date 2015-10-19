function [] = hw2()
clear;
data = load('2015Assign1_1.txt');
p_theta = 0;
case_Q = 'nugget';
t = data(:,1);
y = data(:,2);
% number of unknowns and number of observations
m = 400; n = length(y);
% get the integral matrix
H = shaw(m);
% Q0 and R0
switch case_Q
    case 'nugget'
      Q0 = eye(m);
      X = ones(m,1);
    case 'linear'
      Q0 = -abs(repmat(t,1,n)-repmat(t',n,1));
      X = ones(m,1);
    case 'cubic'
      Q0 = (abs(repmat(t,1,n)-repmat(t',n,1))).^3;
      X = [ones(m,1),t];
end
R0 = 1e-6*eye(n);
% PHI
PHI = H*X;
% Transformation
T = null(PHI')';
z = T*y; nz = length(z);
% the negative loglikelihood function L(theta)
    function res = GetL(theta)
        PSI = H*theta(1)*Q0*H'+theta(2)*R0;
        C = T*PSI*T';
%         log_det_TPSI = trace(logm(C));
        log_det_TPSI = sum(log(eig(C)));
        res = 0.5*log_det_TPSI+0.5*z'*inv(C)*z;
    end
%% Question 1
% %{
n_theta = 10;
switch case_Q
    case 'nugget'
      th1_min = 2; th1_max = 30;
    case 'linear'
      th1_min = 0.1; th1_max = 1;
    case 'cubic'
      th1_min = 0.1; th1_max = 2;
end
theta1 = linspace(th1_min,th1_max,n_theta);
theta2 = linspace(0.5,2,n_theta);
L = zeros(n_theta);
X = zeros(n_theta,1);
Y = zeros(n_theta,1);
count = 0;
for i=1:n_theta
    for j = 1:n_theta
        L(j,i) = GetL([theta1(i),theta2(j)]);
        X(i) = theta2(i);
        Y(j) = theta1(j);
    end
    count = count + 1
end
figure(9);
contour(Y,X,L);
fs = 18;
xlabel('{\theta}_1','fontsize',fs);
ylabel('{\theta}_2','fontsize',fs);
colorbar
saveas(gcf,['hw2_q1_contour_' case_Q],'fig');
saveas(gcf,['hw2_q1_contour_' case_Q],'pdf');
saveas(gcf,['hw2_q1_contour_' case_Q],'png');
figure(8);
surf(Y,X,L,L);
xlabel('{\theta}_1','fontsize',fs);
ylabel('{\theta}_2','fontsize',fs);
zlabel('L({\theta})','fontsize',fs);
set(gca,'zdir','reverse');
saveas(gcf,['hw2_q1_surf_' case_Q],'fig');
saveas(gcf,['hw2_q1_surf_' case_Q],'pdf');
saveas(gcf,['hw2_q1_surf_' case_Q],'png');
%}
%% Question 2
% case 1, p(theta) ~ 1
% search with fminsearch
%{
options = optimset('Display','iter','MaxIter',20,...
    'MaxFunEvals',40);
[X,FVAL] = fminsearch(@GetL,[1;1],options);
%}
% search with fmincon
%{
options = optimoptions('fmincon','display','iter',...
    'MaxFunEvals',40);
[X,FVAL,EXITFLAG] = fmincon(@GetL,[1;1],[],[],[],[],[0;0],[100,100],...
    [],options);
%}
CT = zeros(nz,nz,2);
CT(:,:,1) = T*H*Q0*H'*T';
CT(:,:,2) = T*R0*T';
np = 2;
% search with rml
% %{
[theta,V_theta,L_opt] = rml(CT,[1;1],z,1e-6);
theta
%}
% search with fmincon with full Gauss-Newton information
% %{

    function [fun,grad,Hess] = MAP(theta)
        PSI = H*theta(1)*Q0*H'+theta(2)*R0;
        C = T*PSI*T';
        iC = inv(C);
        %         fun = 0.5*trace(logm(C))+...
        fun = 0.5*sum(log(eig(C)))+...
            0.5*z'*iC*z+p_theta*log(theta(1))+p_theta*log(theta(2));
        grad = zeros(np,1); Hess = zeros(np,np);
        for p=1:np
            grad(p) = 0.5*trace(iC*CT(:,:,p))-...
                0.5*z'*iC*CT(:,:,p)*iC*z + p_theta/theta(p);
            for q=1:np
                Hess(p,q) = 0.5*trace(iC*CT(:,:,p)*iC*CT(:,:,q))+...
                    (p==q)*(-p_theta/(theta(p)^2));
            end
        end
    end
options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'GradObj','on','Hessian','user-supplied','Display','iter',...
    'TolFun',1e-6,'MaxFunEvals',40);
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = ...
    fmincon(@MAP,[1;1],[],[],[],[],[0;0],[100;100],...
    [],options);
X
Cov_theta = inv(HESSIAN)
%}
stop =1 ;

end



function Q=linearQ(n,x,var,l)
theta = var/l;
h=abs(repmat(x,1,n)-repmat(x',n,1));
Q=-theta * h;
end



function Q=cubicQ(n,x)
Q=(abs(repmat(x,1,n)-repmat(x',n,1))).^3;
end



