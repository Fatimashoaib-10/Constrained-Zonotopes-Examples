function Xhat = cz_update(Xbar,YK,v)

Sx1 = 5;
Sy1 = 15;
Sx2 = 15;
Sy2 = 15;

n = 4;

if size(Xbar.Z,2)-1 > 3
    ng_desired = 3;
    nc_desired = 2;
    order = (ng_desired-nc_desired)/n;
    
    Xt = reduce(Xbar,'girard',order,nc_desired);
    
    if ~isempty(Xt)
        Xbar = Xt;
    else
        nc_desired = 1;
        Xbar = reduce(Xbar,'girard',order,nc_desired);
    end
end

cf = Xbar.Z(:,1);
Gf = Xbar.Z(:,2:end);
Af = Xbar.A;
bf = Xbar.b;

cf1 = cf(1,:); cf2 = cf(2,:); cf3 = cf(3,:);
Gf1 = Gf(1,:); Gf2 = Gf(2,:); Gf3 = Gf(3,:);

%% zbh

zbh = YK;

ch = zbh.Z(:,1);
Gh = zbh.Z(:,2:end);

%% f(eta) and its Bounds
eta = sym('eta',[size(Gf,2) 1]);
inv_eta = interval(-ones(size(Gf,2),1),ones(size(Gf,2),1));

func = {@(eta) sqrt((Sx1-(cf1+Gf1*eta)).^2+(Sy1-(cf2+Gf2*eta)).^2) + v(1);
    @(eta) (cf3+Gf3*eta) - atan((Sy1-(cf2+Gf2*eta))/(Sx1-(cf1+Gf1*eta))) + v(2);
    @(eta) sqrt((Sx2-(cf1+Gf1*eta)).^2+(Sy2-(cf2+Gf2*eta)).^2) + v(3);
    @(eta) (cf3+Gf3*eta) - atan((Sy2-(cf2+Gf2*eta))/(Sx2-(cf1+Gf1*eta))) + v(4)};

n = size(func,1);
grad{n} = [];
gradfunc{n} = [];
J{n} = [];
fupp{n} = [];
floww{n} = [];
HH{n} = [];

for i = 1:n
    
    grad{i} = jacobian(func{i}(eta),eta);
    gradfunc{i} = matlabFunction(grad{i},'Vars',{eta});
    
    J{i} = interval(gradfunc{i}(inv_eta));
    
    [~,~,fupp{i},floww{i},HH{i}] = decomp_signstable_modified(func{i},...
        inv_eta.sup,inv_eta.inf,J{i}.sup,J{i}.inf);
    
end

%% Combinations of H matrices
[~,ind1]=min(fupp{1}-floww{1});
[~,ind2]=min(fupp{2}-floww{2});
[~,ind3]=min(fupp{3}-floww{3});
[~,ind4]=min(fupp{4}-floww{4});

Q      = [HH{1}(ind1,:);HH{2}(ind2,:);HH{3}(ind3,:);HH{4}(ind4,:)];
pd_up  = [fupp{1}(ind1,:);fupp{2}(ind2,:);fupp{3}(ind3,:);fupp{4}(ind4,:)];
pd_low = [floww{1}(ind1,:);floww{2}(ind2,:);floww{3}(ind3,:);floww{4}(ind4,:)];

cfu = cf;
Gfu = [Gf zeros(size(Gf,1),size(Gh,2)+size(Gh,1))];

Afu = [Af zeros(size(Af,1),size(Gh,2)+size(Gh,1)); Q -Gh 1/2*diag(pd_up-pd_low)];
bfu = [bf;ch-1/2*(pd_low+pd_up)];

czu = conZonotope(cfu,Gfu,Afu,bfu);

czi = removeRedundancies(mptPolytope(czu));

Xhat = conZonotope(czi);