%% Rego 2019 Example 01 (Theorem 02)

clear;close all;%clc
format long g

%% Initial Sets

c  = [-1;1];   % Center
G  = [0.2 0.4 0.2; 0.2 0 -0.2]; % Generators
A  = [2 2 2];  % Constraints Matrix
b  = -3;       % Constraints Vector
X0 = conZonotope(c,G,A,b);

% figure(1);
% plotZono(X0);

ng = size(G,2);  % Number of Generators
nc = length(b);  % Number of Constraints

W = [0;0];
W = conZonotope(W);

%% Nonlinear Dynamics

func = @(x,w) [3.*x(1) - x(1).^2/7 - 4.*x(1).*x(2)./(4+x(1)) + w(1);... % Eq1
                -2.*x(2) + 3.*x(1).*x(2)./(4+x(1)) + w(2)];             % Eq2

number_of_states = 2;

%% State Estimation

% Prediction
Xest  = prediction(func,X0,W,number_of_states,'C1','J1');
Xest2 = prediction(func,X0,W,number_of_states,'C2','J1');
Xest3 = prediction(func,X0,W,number_of_states,'C1','J2');
Xest4 = prediction(func,X0,W,number_of_states,'C2','J2');

figure(1)
plot(Xest,[1,2],'g','Filled',true);grid;
xlim([-5.5 0.5]);ylim([-6.5 -1]);shg
xlabel('x1');ylabel('x2');

figure(1);hold on
plot(Xest2,[1,2],'r','Filled',true);

figure(2)
plot(Xest3,[1,2],'b','Filled',true);grid;
xlim([-5.5 0.5]);ylim([-6.5 -1]);shg
xlabel('x1');ylabel('x2');

figure(2);hold on
plot(Xest4,[1,2],'m','Filled',true);
