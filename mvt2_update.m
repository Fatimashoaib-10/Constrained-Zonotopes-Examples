function Xhat = mvt2_update(Xbar,Yk,v)

Sx1 = 5;
Sy1 = 15;
Sx2 = 15;
Sy2 = 15;

%% Zbh

zbf = zonoBundle(Xbar);

zbh = zonoBundle(Yk);
ch = zbh.Z{1}.center;
Gh = zbh.Z{1}.generators;

%%
for sz = 1:size(zbf.Z,1)
    
    cf = zbf.Z{sz,1}.Z(:,1);
    Gf = zbf.Z{sz,1}.Z(:,2:end);
    
    cf1 = cf(1,:); cf2 = cf(2,:);cf3 = cf(3,:);
    Gf1 = Gf(1,:); Gf2 = Gf(2,:);Gf3 = Gf(3,:);
    
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
    
    tilde_mu0 = cellfun(@(c) c(center(inv_eta)),func,'UniformOutput',false);
    tilde_mu0 = cell2mat(tilde_mu0)-midJ*center(inv_eta);
    
    cfu = cf;
    Gfu = [Gf zeros(size(Gf,1),size(Gh,2)+size(Gh,1))];
    Afu = [midJ -Gh Gm];
    
    cs = vertices(inv_eta);
    cz{size(cs,2)+1} = [];
    
    for i = 0:size(cs,2)
        if (i>0)
            cvertex = cs(:,i);
            tilde_mu0 = cellfun(@(c) c(cvertex),func,'UniformOutput',false);
            tilde_mu0 = cell2mat(tilde_mu0)-midJ*cvertex;
        end
        bfu = [ch-tilde_mu0];
        cz{i+1}= conZonotope(cfu,Gfu,Afu,bfu);
    end
    
    czi = mptPolytope(cz{1});
    
    for i = 2:length(cz)
        czi = removeRedundancies(and(mptPolytope(cz{i}),czi));
    end
    
    if sz == 1
        poly_all = czi;
    else
        poly_all = and(czi,poly_all);
    end
    
%     zb = zonoBundle(czi);
%     zb = reduce(zb,'girard',1);
%     
%     zbi{sz} = zonoBundle(zb);
%     zbi{sz} = zonoBundle(polytope(zbi{sz}));
%     
%     if sz == 1
%         poly_all = polytope(zbi{sz});
%     else
%         poly_all = and(polytope(zbi{sz}),poly_all);
%     end
    
end
Xhat = poly_all;