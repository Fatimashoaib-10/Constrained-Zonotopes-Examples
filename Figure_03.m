%% Rego 2019 (Figure 3)

% function [] = Figure_03(fid)

%% Dynamics

%%% Comment this out and uncomment 'function []=Figure_03(iii)' if you would like to use run_multiple.m
clear;close all;clc
format long g

fid = 8; % Good seeds: 6, 8, 13, 16, 18

rng(fid); 

number_of_states = 2;

func = @(x,w) [3.*x(1) - x(1).^2/7 - 4.*x(1).*x(2)./(4+x(1)) + w(1);... % Eq1
                -2.*x(2) + 3.*x(1).*x(2)./(4+x(1)) + w(2)];             % Eq2

% Output Equation
C  = [1 0;-1 1];
Du = zeros(2,1);
Dv = eye(2);

Outfunc = @(x,v) C*x + Dv*v;

%% Process and Measurement Disturbances
V = conZonotope(interval([-0.4;-0.4],[0.4;0.4]));
W = conZonotope([0;0]);

%% Initial Conditions

x0 = [0.8;0.65];
w0 = [0;0];
% v0 = [0.2;-0.4];
u0 = 0;

X0 = conZonotope([[0.5;0.5],[0.1 0.2 -0.1;0.1 0.1 0.0]]);

%% 
steps = 5;
v = 0.4*(rand(size(C,1),steps)*2-1);

% Memory Allocation
xk{steps}    = [];
yk{steps}    = [];
Xbar{steps}  = []; % Prediction from Rego Algorithm (h from C2)
Xhat{steps}  = []; % Update from Rego Algorithm
Xbar2{steps} = []; % Prediction from Rego Algorithm (J from decomposition)
Xhat2{steps} = []; % Update from Rego Algorithm
Xbar3{steps} = []; % Prediction from Decomposition with Zonotope Bundles
Xhat3{steps} = []; % Update from ZB
Xbar4{steps} = []; % Prediction from Decomposition with Constrained Zonotopes
Xhat4{steps} = []; % Update from CZ
Xbar5{steps} = []; % Combination
Xhat5{steps} = []; % Combination
xs{steps}    = [];

% Initialization
xk{1}    = x0;
yk{1}    = Outfunc(xk{1},v(:,1));
YK       = (yk{1}-Du*u0)+(-Dv*V);

Xbar{1}  = X0;    
Xhat{1}  = update_state(Xbar{1},yk{1},C,Du,u0,Dv,V);
Xbar2{1} = Xbar{1};
Xhat2{1} = Xhat{1};
Xbar3{1} = Xbar{1};
Xhat3{1} = zb_decomp_update(Xbar3{1},YK,v(:,1)); 
Xbar4{1} = Xbar{1};
Xhat4{1} = cz_decomp_update(Xbar4{1},YK,v(:,1)); 
Xbar5{1} = Xbar{1};
Xhat5{1} = Xhat{1};

I0 = zonotope(interval(Xhat{1}));

rng(10); 
size_samples = 100; % Increasing samples will take a long time
xi = (-1+2.*rand(2,size_samples));
xs0 = [];

for i=1:length(xi)
    x=I0.Z*[1;xi(:,i)];
    if in(Xhat{1},x)
        xs0=[xs0,  x];
    end
end

xs{1} = xs0;

figure(fid);hold on;grid off;box on;
plot(Xhat{1},[1,2],'g');%,'Filled',true);
plot(Xhat2{1},[1,2],'r');%,'Filled',true);
plot(Xhat3{1},[1,2],'b');%,'Filled',true);
plot(Xhat4{1},[1,2],'m');%,'Filled',true);
plot(Xhat5{1},[1,2],'c');%,'Filled',true);
scatter(xs{1}(1,:),xs{1}(2,:),10,'ko')
drawnow;
%%
for k = 2:steps
    xk{k} = func(xk{k-1},w0);
    yk{k} = Outfunc(xk{k},v(:,k));
    YK = (yk{k}-Du*u0)+(-Dv*V);
    
    xs{k} = zeros(number_of_states ,0);
    
    for ii=1:size(xs{k-1},2)
        xt = func(xs{k-1}(:,ii),w0);
        if in(YK,C*xt)
            xs{k} = [xs{k},xt];%func(xs{k-1}(:,ii),w0)];
        end
    end
    
    % Rego Theorem 02 (Mean-Value Extension)
    tic
    Xbar{k} = prediction(func,Xhat{k-1},W,number_of_states,'C2','J1');
    Xhat{k} = update_state(Xbar{k},yk{k},C,Du,u0,Dv,V);
    toc
    
    % Rego Theorem 02 (Mean-Value Extension with Decomposition-based J)
    tic
    Xbar2{k} = prediction(func,Xhat2{k-1},W,number_of_states,'C2','J2');
    Xhat2{k} = update_state(Xbar2{k},yk{k},C,Du,u0,Dv,V);
    toc
    
    % Decomposition with Zonotope Bundles
    tic
    Xbar3{k} = zb_decomp_prediction(Xhat3{k-1});
    Xhat3{k} = zb_decomp_update(Xbar3{k},YK,v(:,k));
    toc
    
    % Decomposition with Constrained Zonotopes
    tic
    Xbar4{k} = cz_decomp_prediction(Xhat4{k-1});
    Xhat4{k} = cz_decomp_update(Xbar4{k},YK,v(:,k));
    toc
    
    % Combined
    tic
    Xbar5{k} = mptPolytope(prediction(func,Xhat5{k-1},W,number_of_states,'C2','J1'));
    Xbar5{k} = removeRedundancies(and(Xbar5{k},mptPolytope(prediction(func,Xhat5{k-1},W,number_of_states,'C2','J2'))));
    Xbar5{k} = removeRedundancies(and(Xbar5{k},mptPolytope(zb_decomp_prediction(Xhat5{k-1}))));
    Xbar5{k} = removeRedundancies(and(Xbar5{k},mptPolytope(cz_decomp_prediction(Xhat5{k-1}))));
    Xbar5{k} = conZonotope(Xbar5{k});
    Xhat5{k} = update_state(Xbar5{k},yk{k},C,Du,u0,Dv,V);
    toc
    
    % Plots
    plot(Xhat{k},[1,2],'g');
    plot(Xhat2{k},[1,2],'r');
    plot(Xhat3{k},[1,2],'b');
    plot(Xhat4{k},[1,2],'m');
    plot(Xhat5{k},[1,2],'c');
    scatter(xs{k}(1,:),xs{k}(2,:),15,'ko');
    drawnow;
end

title('Figure 3 (Rego)')
xlabel('x1')
ylabel('x2')
legend('Rego','J from decomposition','decomposition with ZB','decomposition with CZ','Combined')
legend boxoff