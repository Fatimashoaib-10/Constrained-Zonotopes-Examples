function Xbar = zb_decomp_prediction(X0)

Z0 = zonoBundle(X0);

for sz = 1:size(Z0.Z,1)
    
    c = Z0.Z{sz,1}.Z(:,1); c1 = c(1,:); c2 =c(2,:);
    G = Z0.Z{sz,1}.Z(:,2:end); G1 = G(1,:); G2 =G(2,:);
    
    %% f(eta) and its Bounds
    
    n   = 2; % number of states
    eta = sym('eta',[size(G,2) 1]);
    inv_eta = interval(-ones(size(G,2),1),ones(size(G,2),1));
    
    func = {@(eta)  3.*(c1+G1*eta) - (c1+G1*eta).^2/7 - 4.*(c1+G1*eta).*(c2+G2*eta)./(4+(c1+G1*eta));...
        @(eta) -2.*(c2+G2*eta) + 3.*(c1+G1*eta).*(c2+G2*eta)./(4+(c1+G1*eta))};
    
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
    H         = [HH{1}(ind1,:);HH{2}(ind2,:)];
    gd_up     = [fupp{1}(ind1,:);fupp{2}(ind2,:)];
    gd_low    = [floww{1}(ind1,:);floww{2}(ind2,:)];
    
    %% Zonotope Bundle

    cf = 1/2*(gd_low + gd_up);
    Gf = [H 1/2*diag(gd_up - gd_low)];
    
    zf = zonotope(cf,Gf);
    zb   = zonoBundle(zf);
    
    zbi{sz} = zonoBundle(zb);
    zbi{sz} = zonoBundle(polytope(zbi{sz}));
    
    if sz == 1
        poly_all = polytope(zbi{sz});
    else
        poly_all = and(polytope(zbi{sz}),poly_all);
    end
end

Xbar = poly_all;