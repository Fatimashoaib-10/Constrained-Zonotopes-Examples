function Xbar = cz_prediction(X0,w)

n = 3;
if size(X0.Z,2)-1 > 3
    ng_desired = 3;
    nc_desired = 2;
    order = (ng_desired-nc_desired)/n;
    
    Xt = reduce(X0,'girard',order,nc_desired);
    
    if ~isempty(Xt)
        X0 = Xt;
    else
        nc_desired = 1;
        X0 = reduce(X0,'girard',order,nc_desired);
    end
end

c = X0.Z(:,1);
G = X0.Z(:,2:end);
A = X0.A;
b = X0.b;

c1 = c(1,:); c2 = c(2,:); c3 = c(3,:);
G1 = G(1,:); G2 = G(2,:); G3 = G(3,:);

%% f(eta) and its Bounds

eta = sym('eta',[size(G,2) 1]);
inv_eta = interval(-ones(size(G,2),1),ones(size(G,2),1));

phi_omega = 0.30; % Displacement Velocity
phi_theta = 0.15; % Angular Velocity
T0        = 1.00; % Sampling Period

func = {@(eta)  (c1+G1*eta) + T0.*phi_omega.*cos((c3+G3*eta)) + w(1);...
    @(eta)  (c2+G2*eta) + T0.*phi_omega.*sin((c3+G3*eta)) + w(2);...
    @(eta)  (c3+G3*eta) + T0.*phi_theta + w(3)};

n   = size(func,1);

grad{n}     = [];
gradfunc{n} = [];
J{n}        = [];
fupp{n}     = [];
floww{n}    = [];
HH{n}       = [];

for i = 1:n
    
    grad{i}      = jacobian(func{i}(eta),eta);
    gradfunc{i}  = matlabFunction(grad{i},'Vars',{eta});
    
    J{i}  = interval(gradfunc{i}(inv_eta));
    
    [~,~,fupp{i},floww{i},HH{i}] = decomp_signstable_modified(func{i},...
        inv_eta.sup,inv_eta.inf,J{i}.sup,J{i}.inf);
end

%% Combinations of H matrices
[~,ind1]  = min(fupp{1}-floww{1});
[~,ind2]  = min(fupp{2}-floww{2});
[~,ind3]  = min(fupp{3}-floww{3});

H      = [HH{1}(ind1,:);HH{2}(ind2,:);HH{3}(ind3,:);];
gd_up  = [fupp{1}(ind1,:);fupp{2}(ind2,:);fupp{3}(ind3,:)];
gd_low = [floww{1}(ind1,:);floww{2}(ind2,:);floww{3}(ind3,:)];

%% Constrained Zonotopes
cf = 1/2*(gd_low + gd_up);
Gf = [H 1/2*diag(gd_up - gd_low)];
Af = [A zeros(size(A,1),size(Gf,2)-size(A,2))];
bf = b;

cz = conZonotope(cf,Gf,Af,bf);

%% Intersection

czi = removeRedundancies(mptPolytope(cz));

Xbar = conZonotope(czi);