function Xbar = mvt_prediction(X0,w)

n = 3;

% X0 = conZonotope(mptPolytope(X0));

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


for i = 1:n
    
    grad{i}      = jacobian(func{i}(eta),eta);
    gradfunc{i}  = matlabFunction(grad{i},'Vars',{eta});
    J{i}  = interval(gradfunc{i}(inv_eta));
    
    for j = 1:size(grad{i},2)
        f{i,j}   = matlabFunction(grad{i}(1,j),'Vars',{eta});
        df{i,j}  = matlabFunction(jacobian(grad{i}(1,j),eta),'Vars',{eta});
        dfb{i,j} = interval(df{i,j}(inv_eta));
        [Fmax(i,j),Fmin(i,j)] = decomp_signstable(f{i,j},...
            inv_eta.sup,inv_eta.inf,...
            dfb{i,j}.sup,dfb{i,j}.inf);
        J2(i,j) = interval(Fmin(i,j),Fmax(i,j));
    end
    
end

J = J2;

midJ = center(J);

Gm = zeros(n,1);

for i = 1:size(J,1)
    for j = 1:size(J,2)
        Gm(i) = Gm(i)+rad(J(i,j));
    end
end

Gm = diag(Gm);

tilde_mu0 = cellfun(@(c) c(center(inv_eta)),func,'UniformOutput',false);
cf = cell2mat(tilde_mu0)-midJ*center(inv_eta);

cs = vertices(inv_eta);
cz{size(cs,2)} = [];

Gf = [midJ Gm];
Af = [A zeros(size(A,1),size(Gm,2))];

for i = 0:size(cs,2)
    if (i>0)
        cvertex = cs(:,i);
        tilde_mu = cellfun(@(c) c(cvertex),func,'UniformOutput',false);
        cf = cell2mat(tilde_mu)-midJ*cvertex;
    end
    cz{i+1} = conZonotope(cf,Gf,Af,b);
    
end

czi = mptPolytope(cz{1});

for i = 2:length(cz)
    czi = removeRedundancies(and(mptPolytope(cz{i}),czi));
end

Xbar = conZonotope(czi);