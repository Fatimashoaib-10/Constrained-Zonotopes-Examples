function Xhat = mvx_update(Xbar,YK,v)

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

tilde_mu0 = cellfun(@(c) c(zeros(size(eta,1),1)),func,'UniformOutput',false);
tilde_mu0 = cell2mat(tilde_mu0);

cfu = cf;
Gfu = [Gf zeros(size(Gf,1),size(Gh,2)+size(Gh,1))];
Afu = [Af zeros(size(Af,1),size(Gh,2)+size(Gh,1)); midJ -Gh Gm];
bfu = [bf;ch-tilde_mu0];

Xhat = conZonotope(cfu,Gfu,Afu,bfu);
