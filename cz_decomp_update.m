function Xhat = cz_decomp_update(Xbar,YK,v)

n = 2; % number of states

cf = Xbar.Z(:,1);
Gf = Xbar.Z(:,2:end);
Af = Xbar.A;
bf = Xbar.b;

cf1 = cf(1,:); cf2 = cf(2,:);
Gf1 = Gf(1,:); Gf2 = Gf(2,:);

%% zbh
zbh = YK;

ch = zbh.Z(:,1);
Gh = zbh.Z(:,2:end);

%% f(eta) and its Bounds
eta = sym('eta',[size(Gf,2) 1]);
inv_eta = interval(-ones(size(Gf,2),1),ones(size(Gf,2),1));

func = {@(eta) cf1+Gf1*eta + v(1);
    @(eta) -(cf1+Gf1*eta) + (cf2+Gf2*eta) + v(2)};

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

[~,ind1] = min(fupp{1}-floww{1});
[~,ind2] = min(fupp{2}-floww{2});

Q        = [HH{1}(ind1,:);HH{2}(ind2,:)];
pd_up    = [fupp{1}(ind1,:);fupp{2}(ind2,:)];
pd_low   = [floww{1}(ind1,:);floww{2}(ind2,:)];

cfu = cf;
Gfu = [Gf zeros(size(Gf,1),size(Gh,2)+size(Gh,1))];
Afu = [Af zeros(size(Af,1),size(Gh,2)+size(Gh,1));Q -Gh 1/2*diag(pd_up-pd_low)];
bfu = [bf;ch-1/2*(pd_low+pd_up)];

czu = conZonotope(cfu,Gfu,Afu,bfu);

Xhat = czu;