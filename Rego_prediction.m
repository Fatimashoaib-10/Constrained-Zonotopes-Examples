function X_predicted = Rego_prediction(func,X0,w,W,number_of_states)

% Theorem 02

% if size(X0.Z,2)-1 > 3
%     ng_desired = 3;
%     nc_desired = 2;
%     order = (ng_desired-nc_desired)/3;
%     
%     Xt = reduce(X0,'girard',order,nc_desired);
%     
%     if ~isempty(Xt)
%         X0 = Xt;
%     else
%         nc_desired = 1;
%         X0 = reduce(X0,'girard',order,nc_desired);
%     end
% end

% Gradient
x = sym('x',[number_of_states 1]);
grad     = jacobian(matlabFunction(func(x,w)),x);
gradfunc = matlabFunction(grad,'Vars',{x});    % Gradient

% Interval Hulls of X and W

interval_X = interval(X0);
interval_W = interval(W);

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

Z = func(h,center(interval_W));
Z = conZonotope(Z);

% Predicted Set
X_predicted = Z + theorem_01(J,X0+(-h));

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