%% Rego 2019 Example 01 (Update)

clear;close all;clc
format short g

%% Initial Conditions
X0 = conZonotope([[0.5;0.5],[0.1 0.2 -0.1;0.1 0.1 0.0]]);
V0 = conZonotope(interval([0;0],[0.4;0.4]));
x0 = [0.8;0.65];
v0 = [0.17;0.02];
u0 = 0;

% yk = ck*xk + Dv*vk
C  = [1 0;-1 1];
Du = zeros(2,1);
Dv = eye(2);

xk = x0;
vk = v0;
uk = u0;

yk   = C*xk + Dv*vk; % Measurement
Xbar = X0; % X_predicted
V    = V0;

Xhat = update_state(Xbar,yk,C,Du,uk,Dv,V); %X_updated

figure(1)
plot(Xbar,[1,2],'r','Filled',true);grid
xlim([0 1]);ylim([0 1]);
xlabel('x1');ylabel('x2');shg
hold on

figure(1)
plot(Xhat,[1,2],'g','Filled',true);
lgd = legend('$\bar{X}$','$\hat{X}$');
set(lgd,'Interpreter', 'latex');
lgd.FontSize = 14;


