function X_updated = Rego_update(func,X0,YK,v,n)

% Theorem 02
cf = X0.Z(:,1);
Gf = X0.Z(:,2:end);
Af = X0.A;
bf = X0.b;

ch = YK.Z(:,1);
Gh = YK.Z(:,2:end);

% Gradient
x = sym('x',[n 1]);
grad     = jacobian(matlabFunction(func(x,v)),x);

% Interval Hulls of X and W

interval_X = interval(X0);

mJ        = size(grad,1);
nJ        = size(grad,2);
f         = cell(mJ,nJ);
df        = cell(mJ,nJ);
dfb       = cell(mJ,nJ);
Fmin      = zeros(mJ,nJ);
Fmax      = zeros(mJ,nJ);
J(mJ,nJ) = interval(0,0);

for i = 1:mJ
    for j = 1:nJ
        f{i,j}   = matlabFunction(grad(i,j),'Vars',{x});
        df{i,j}  = matlabFunction(jacobian(grad(i,j),x),'Vars',{x});
        dfb{i,j} = interval(df{i,j}(interval_X));
        [Fmax(i,j),Fmin(i,j)] = decomp_signstable(f{i,j},...
            interval_X.sup,interval_X.inf,...
            dfb{i,j}.sup,dfb{i,j}.inf);
        J(i,j) = interval(Fmin(i,j)-eps,Fmax(i,j)+eps);
    end
end

h = center(interval_X);

Z = func(h,v);
Z = conZonotope(Z);

% Predicted Set
Xh = Z + theorem_01(J,X0+(-h));

cfu = cf;
Gfu = [Gf zeros(size(Gf,1),size(Gh,2)+size(Gh,1))];
Afu = [Af zeros(size(Af,1),size(Gh,2)+size(Gh,1)); Xh.Z(:,2:end) -Gh];
bfu = [bf;ch-Xh.Z(:,1)];

X_updated = conZonotope(cfu,Gfu,Afu,bfu);

end

function S = theorem_01(J,X)

% Constraint Reduction (X0 -> Xbar)
Xbar = reduce(X,'girard',1,0);

pbar = Xbar.Z(:,1);     % Center of Reduced Zonotope
Mbar = Xbar.Z(:,2:end); % Generators of Reduced Zonotope

% P Matrix
mgbar = size(Mbar,1);
ngbar = size(Mbar,2);
n = size(J,1);

m = (J-center(J))*pbar;

P = zeros(1,n);
for i = 1:n
    term = 0;
    for j = 1:ngbar
        for k = 1:mgbar
            term = term+((2*rad(J(i,k)*abs(Mbar(k,j)))));
        end
    end
    P(i) = rad(m(i)) + 1/2*term;
end

P = diag(P);

% S (Output of Theorem 01) -> S = center(J)*X + P B_inf^n
S1  = center(J)*X;

Sc = S1.Z(:,1);
SG = [S1.Z(:,2:end) P];

SA = [S1.A zeros(size(S1.A,1),n)];
Sb = S1.b;
S = conZonotope(Sc,SG,SA,Sb);

end