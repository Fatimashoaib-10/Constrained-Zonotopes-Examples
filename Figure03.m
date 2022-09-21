%% Unicycle

clear;close all;%clc
clear functions
format short g
rng(3,'twister');
warning('off');
set(groot,'defaultlineLineWidth',2)

steps = 5;

flagZB       = 1;
flagCZ       = 0;
flagMVT      = 0;
flagREGO     = 1;
flagREGOMVT  = 0;
flagCombined = 0;

% MVT_prediction = @mvx_prediction;
% MVT_update     = @mvx_update;
MVT_prediction = @mvt2_prediction;
MVT_update     = @mvt2_update;

%% Nonlinear Dynamics

number_of_states = 3;

phi_omega = 0.30; % Displacement Velocity
phi_theta = 0.15; % Angular Velocity
T0        = 1.00; % Sampling Period

func = @(x,w) [x(1) + T0.*phi_omega.*cos(x(3)) + w(1);...
    x(2) + T0.*phi_omega.*sin(x(3)) + w(2);...
    x(3) + T0.*phi_theta + w(3)];

% Ouput Equation
Sx1 = 5;
Sy1 = 15;
Sx2 = 15;
Sy2 = 15;

Dv = eye(4);

Outfunc = @(x,v) [sqrt((Sx1-x(1)).^2+(Sy1-x(2)).^2)+v(1);
    x(3) - atan((Sy1-x(2))/(Sx1-x(1))) + v(2);
    sqrt((Sx2-x(1)).^2+(Sy2-x(2)).^2)+v(3);
    x(3) - atan((Sy2-x(2))/(Sx2-x(1))) + v(4)];

%% Process & Measurement Disturbances

w0 = [0.5*rand(1,steps)-0.3;
    0.3*rand(1,steps)-0.2;
    0.6*rand(1,steps)-0.4];

w0 = 0.2*eye(number_of_states)*w0;

W = conZonotope(interval([0.5*0-0.3;0.3*0-0.2;0.6*0-0.4],...
    [0.5*1-0.3;0.3*1-0.2;0.6*1-0.4]));

v0 = [0.02*rand(1,steps)-0.01;
    0.03*rand(1,steps)-0.01;
    0.03*rand(1,steps)-0.02;
    0.05*rand(1,steps)-0.03];

V = conZonotope(interval([0.02*0-0.01;0.03*0-0.01;0.03*0-0.02;0.05*0-0.03],...
    [0.02*1-0.01;0.03*1-0.01;0.03*1-0.02;0.05*1-0.03]));

WM = conZonotope(interval([0;0;0],[0;0;0]));

v0 = v0*0;
w0 = w0*0;

%% Initial Conditions

x0 = [0.1;0.2;1];

percentage = [10;10;20];

x01 = [x0(1)*(1-sign(x0(1))*percentage(1)/100),x0(1)*(1+sign(x0(1))*percentage(1)/100)];
x02 = [x0(2)*(1-sign(x0(2))*percentage(2)/100),x0(2)*(1+sign(x0(2))*percentage(2)/100)];
x03 = [x0(3)*(1-sign(x0(3))*percentage(3)/100),x0(3)*(1+sign(x0(3))*percentage(3)/100)];

X0 = interval([x01(1);x02(1);x03(1)],[x01(2);x02(2);x03(2)]);
X0 = conZonotope(X0);

%% State Estimation

xk{steps}           = [];
yk{steps}           = [];
YK{steps}           = [];

if(flagZB)
    XbarZB{steps}       = [];
    XhatZB{steps}       = [];
end

if(flagCZ)
    XbarCZ{steps}       = [];
    XhatCZ{steps}       = [];
end

if(flagMVT)
    XbarMVT{steps}      = [];
    XhatMVT{steps}      = [];
end

if(flagREGO)
    XbarREGO{steps}     = [];
    XhatREGO{steps}     = [];
end

if(flagREGOMVT)
    XbarREGOMVT{steps}  = [];
    XhatREGOMVT{steps}  = [];
end

if(flagCombined)
    XbarCombined{steps} = [];
    XhatCombined{steps} = [];
end

xs{steps}           = [];

xk{1} = x0;
yk{1} = Outfunc(xk{1},v0(:,1));
YK{1} = yk{1}+(-Dv*V);

if(flagZB)
    fprintf('\n ZB: %d\n',1);
    XbarZB{1}       = X0;
    tic;XhatZB{1}   = zb_update(XbarZB{1},YK{1},v0(:,1));toc
end

if(flagCZ)
    fprintf('\n CZ: %d\n',1);
    XbarCZ{1}       = X0;
    tic;XhatCZ{1}   = cz_update(XbarCZ{1},YK{1},v0(:,1));toc
end

if(flagMVT)
    fprintf('\n MVT: %d\n',1);
    XbarMVT{1}      = X0;
    tic;XhatMVT{1}  = MVT_update(XbarMVT{1},YK{1},v0(:,1));toc
end

if(flagREGO)
    fprintf('\n REGO: %d\n',1);
    XbarREGO{1}     = X0;
    tic;XhatREGO{1} = Rego_update(Outfunc,XbarREGO{1},YK{1},v0(:,1),number_of_states);toc;
end

if(flagREGOMVT)
    fprintf('\n REGO+MVT: %d\n',1);
    XbarREGOMVT{1}  = X0;
    if(flagMVT)
        XhatREGOMVT{1}  = XhatMVT{1};
    else
        XhatREGOMVT{1}  = MVT_update(XbarREGOMVT{1},YK{1},v0(:,1));
    end
end
clear functions
if(flagCombined)
    fprintf('\n Combined: %d\n',1);
    XbarCombined{1} = X0;
    
    if(flagZB);c{1} = XhatZB{1};end
    if(flagCZ);c{2} = XhatCZ{1};end
    if(flagMVT);c{3} = XhatMVT{1};end
    if(flagREGO);c{4} = XhatREGO{1};end
    
    c = c(~cellfun('isempty',c));
    
    XhatCombined{1} = mptPolytope(c{1});
    for m = 2:length(c)
        XhatCombined{1} = removeRedundancies(and(mptPolytope(c{m}),XhatCombined{1}));
    end
    XhatCombined{1} = conZonotope(XhatCombined{1});
end
clear functions
I0 = zonotope(interval(XhatZB{1}));

size_samples = 100; % Increasing samples will take a long time
xi = (-1+2.*rand(number_of_states,size_samples));
xs0 = [];

for i=1:length(xi)
    x=I0.Z*[1;xi(:,i)];
    if in(XhatZB{1},x)
        xs0 = [xs0,  x];
    end
end

xs{1} = xs0;

figure(1);hold on;grid on
if(flagZB);plot(XhatZB{1},[1,2],'r');end
if(flagCZ);plot(XhatCZ{1},[1,2],'g');end
if(flagMVT);plot(XhatMVT{1},[1,2],'b');end
if(flagREGO);plot(XhatREGO{1},[1,2],'color',[0.9290, 0.6940, 0.1250]);end
if(flagREGOMVT);plot(XhatREGOMVT{1},[1,2],'c');end
if(flagCombined);plot(XhatCombined{1},[1,2],'m');end

scatter(xs{1}(1,:),xs{1}(2,:),10,'ko');
drawnow;

%%
for k = 2:steps
    
    xk{k} = func(xk{k-1},w0(:,k));
    yk{k} = Outfunc(xk{k},v0(:,k));
    YK{k} = yk{k}+(-Dv*V);
    
    xs{k} = zeros(number_of_states ,0);
    
    for ii=1:size(xs{k-1},2)
        xt=func(xs{k-1}(:,ii),w0(:,k));
        if in(YK{k},Outfunc(xt,v0(:,k)))
            xs{k} = [xs{k},xt];
        end
    end
    
    if(flagZB)
        % Zonotope Bundle Method
        fprintf('\n ZB: %d\n',k);
        tic;XbarZB{k} = zb_prediction(XhatZB{k-1},w0(:,k));toc;
        tic;XhatZB{k} = zb_update(XbarZB{k},YK{k},v0(:,k));toc;
    end
    clear functions
    if(flagCZ)
        % Constrained Zonotopes Method
        fprintf('\n CZ: %d\n',k);
        tic;XbarCZ{k} = cz_prediction(XhatCZ{k-1},w0(:,k));toc
        tic;XhatCZ{k} = cz_update(XbarCZ{k},YK{k},v0(:,k));toc;
    end
    clear functions
    if(flagMVT)
        % Mean Value Theorem
        fprintf('\n MVT: %d\n',k);
        XbarMVT{k} = MVT_prediction(XhatMVT{k-1},w0(:,k));toc
        XhatMVT{k} = MVT_update(XbarMVT{k},YK{k},v0(:,k));toc;
    end
    clear functions
    if(flagREGO)
        % Rego Prediction Update
        fprintf('\n REGO: %d\n',k);
        tic;XhatREGO{k-1}   = conZonotope(mptPolytope(XhatREGO{k-1}));
        XbarREGO{k} = Rego_prediction(func,XhatREGO{k-1},w0(:,k),WM,3);toc
        tic;XbarREGO{k}     = conZonotope(mptPolytope(XbarREGO{k}));
        XhatREGO{k} = Rego_update(Outfunc,XbarREGO{k},YK{k},v0(:,k),number_of_states);toc;
    end
    clear functions
    if(flagREGOMVT)
        % Rego Prediction + MVT Update
        fprintf('\n REGO+MVT: %d\n',k);
%         XhatREGOMVT{k-1} = conZonotope(XhatREGOMVT{k-1});
        tic;XbarREGOMVT{k} = Rego_prediction(func,XhatREGOMVT{k-1},w0(:,k),WM,3);toc
        tic;XbarREGOMVT{k} = conZonotope(mptPolytope(XbarREGOMVT{k}));
        XhatREGOMVT{k} = MVT_update(XbarREGOMVT{k},YK{k},v0(:,k));toc;
        XhatREGOMVT{k} = conZonotope(XhatREGOMVT{k});
    end
    clear functions
    if(flagCombined)
        % Combined
        fprintf('\n Combined: %d\n',k);
        tic;
        XhatCombined{k-1} = conZonotope(mptPolytope(XhatCombined{k-1}));
        if(flagZB);a{1}   = zb_prediction(XhatCombined{k-1},w0(:,k));end
        if(flagCZ);a{2}   = cz_prediction(XhatCombined{k-1},w0(:,k));end
        if(flagMVT);a{3}  = MVT_prediction(XhatCombined{k-1},w0(:,k));end
        if(flagREGO)
            a{4} = Rego_prediction(func,XhatCombined{k-1},w0(:,k),WM,3);
            a{4} = conZonotope(mptPolytope(a{4}));
        end
        
        %     a{3} = zonoBundle(mvx_prediction(XhatCombined{k-1},w0(:,k)));
        %     a{4} = zonoBundle(Rego_prediction(func,XhatCombined{k-1},w0(:,k),WM,3));
        
        a = a(~cellfun('isempty',a));
        
        XbarCombined{k} = mptPolytope(a{1});
        for m = 2:length(a)
            XbarCombined{k} = removeRedundancies(and(mptPolytope(a{m}),XbarCombined{k}));
        end
        XbarCombined{k} = conZonotope(XbarCombined{k});
        toc;
        
        tic;
        if(flagZB);b{1}   = zb_update(XbarCombined{k},YK{k},v0(:,k));end
        if(flagCZ);b{2}   = cz_update(XbarCombined{k},YK{k},v0(:,k));end
        if(flagMVT);b{3}  = MVT_update(XbarCombined{k},YK{k},v0(:,k));end
        if(flagREGO)
            XbarCombined{k} = conZonotope(mptPolytope(XbarCombined{k}));
            b{4} = Rego_update(Outfunc,XbarCombined{k},YK{k},v0(:,k),number_of_states);
            b{4} = conZonotope(mptPolytope(b{4}));
        end
       
        b = b(~cellfun('isempty',b));
        
        XhatCombined{k} = mptPolytope(b{1});
        for m = 2:length(b)
            XhatCombined{k} = removeRedundancies(and(mptPolytope(b{m}),XhatCombined{k}));
        end
        XhatCombined{k} = conZonotope(XhatCombined{k});
        
        toc;
    end
    clear functions
    % Plots
    if(flagZB);plot(XhatZB{k},[1,2],'r');end
    if(flagCZ);plot(XhatCZ{k},[1,2],'g');end
    if(flagMVT);plot(XhatMVT{k},[1,2],'b');end
    if(flagREGO);plot(XhatREGO{k},[1,2],'color',[0.9290, 0.6940, 0.1250]);end
    if(flagREGOMVT);plot(XhatREGOMVT{k},[1,2],'c');end
    if(flagCombined);plot(XhatCombined{k},[1,2],'m');end
    scatter(xs{k}(1,:),xs{k}(2,:),10,'ko');
    drawnow;
    
end

xlabel('Sx');ylabel('Sy');title('Unicycle Robot Trajectory');
%%
lgdE = [flagZB;flagCZ;flagMVT;flagREGO;flagREGOMVT;flagCombined];
lgdM = {'Zonotope Bundle Method',...
    'Constrained Zonotope Method',...
    'Mean Value Theorem',...
    'REGO',...
    'REGO Prediction + MVT Update',...
    'Combination'};

for x = 1:length(lgdE)
    if(lgdE(x) == 1)
        Legend{x} = lgdM{x};
    end
end

Legend = Legend(~cellfun('isempty',Legend));
legend(Legend,'Location','best');
