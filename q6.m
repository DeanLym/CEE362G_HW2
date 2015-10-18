function [] = q6(theta1,theta2,model)
close all;
data = load('2015Assign1_1.txt');
t = data(:,1);
y = data(:,2);

n=400;
ind = 20;
l_s = pi*5;
H = shaw(n); 

%% case 1 Nugget Effect
R = theta2*eye(n);
Q = theta1*eye(n);
X = ones(n,1);
[s1,V1,lamda1,mu1] = GenLinInv(y,H,R,X,Q);
% plotresult(t,s,s_true,61);
%% case 2 Linear Covariance
R = theta2*eye(n);
Q = linearQ(n,t,theta1,l_s);
X = ones(n,1);
[s2,V2,lamda2,mu2] = GenLinInv(y,H,R,X,Q);
%% case 3 NonLinear Covariance
R = theta2*eye(n);
Q = cubicQ(n,t,theta1,l_s);
X = zeros(n,2);
X(:,1)=ones(n,1);
X(:,2)=t;
[s3,V3,lamda3,mu3] = GenLinInv(y,H,R,X,Q);

figure(61+ind);
hold on;
if (model==1)
    plot(t,s1,'r','linewidth',2);
    plot(t,s1+2*sqrt(diag(V1)),'r--','linewidth',2);
    plot(t,s1-2*sqrt(diag(V1)),'r--','linewidth',2);

else if (model==2)
        plot(t,s2,'r','linewidth',2);
        plot(t,s1+2*sqrt(diag(V2)),'r--','linewidth',2);
        plot(t,s1-2*sqrt(diag(V2)),'r--','linewidth',2);

    else if (model==3)
            plot(t,s3,'r','linewidth',2);
            plot(t,s1+2*sqrt(diag(V3)),'r--','linewidth',2);
            plot(t,s1-2*sqrt(diag(V3)),'r--','linewidth',2);

        end
    end
end

plot(t,y,'--','linewidth',2);
legend('$\hat{s}$(t) Nugget Covariance','credibility intervals','credibility intervals','y(t)',...
    'Location','NorthEastOutside');
set(legend,'Interpreter','latex');
set(legend,'Fontsize',12);
set(legend,'Fontname','Times New Roman');

%% Question 7
nr=10;
s_c = zeros(n,nr);
figure(65+ind);
% title('Conditional realizations of $\hat{s}$(t)','Interpreter','latex','Fontsize',12,'Fontname','Times New Roman');
hold on;
% Q = theta1*eye(n);
if (model==1)
    Q = theta1*eye(n);
    lamda = lamda1;
else if (model==2)
        Q = linearQ(n,t,theta1,l_s);
        lamda = lamda2;
    else if (model==3)
            Q = cubicQ(n,t,theta1,l_s);
            lamda = lamda3;
        end
    end
end
for i=1:nr
v_u = simv(R);
s_u = simvm2(1,Q,X);
s_c(:,i) = s_u + lamda*(y+v_u-H*s_u);
subplot(2,5,i);
hold on;
plot(t,s_c(:,i));
plot(t,s2,'r');
ylim([-3,3]);
end

end


function [] = plotresult(s)
% ind figure number
figure(ind);
hold on;
plot(t,s,'bo');
plot(t,s_true,'r-','linewidth',2);

end




function Q=linearQ(n,x,var,l)

theta = var/l;
h=abs(repmat(x,1,n)-repmat(x',n,1));
Q=-theta * h;

end

function Q=cubicQ(n,x,var,l)

theta = var/l;
h=abs(repmat(x,1,n)-repmat(x',n,1));
Q=-theta * h^3;

end


function Q=smoothCov(n,t)

h = abs(repmat(t,1,n) - repmat(t',n,1));
Q = h.^3;

end
