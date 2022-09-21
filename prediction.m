function X_predicted = prediction(func,X0,W,number_of_states,h_str,J_str)

% Theorem 02

% Gradient
x = sym('x',[number_of_states 1]);
grad     = jacobian(matlabFunction(func(x,zeros(1,number_of_states))),x);
gradfunc = matlabFunction(grad,'Vars',{x});    % Gradient

% Interval Hulls of X and W

interval_X = interval(X0);
interval_W = interval(W);

J1 = gradfunc(interval_X);

switch (J_str)
    case 'J1'
        J = J1;
    case 'J2'
        mJ        = size(grad,1);
        nJ        = size(grad,2);
        f         = cell(mJ,nJ);
        df        = cell(mJ,nJ);
        dfb       = cell(mJ,nJ);
        Fmin      = zeros(mJ,nJ);
        Fmax      = zeros(mJ,nJ);
        J2(mJ,nJ) = interval(0,0);
        
        for i = 1:mJ
            for j = 1:nJ
                f{i,j} = matlabFunction(grad(i,j),'Vars',{x});
                df{i,j} = matlabFunction(jacobian(grad(i,j),x),'Vars',{x});
                dfb{i,j} = df{i,j}(interval_X);
                [Fmax(i,j),Fmin(i,j)] = decomp_signstable(f{i,j},...
                    interval_X.sup,interval_X.inf,...
                    dfb{i,j}.sup,dfb{i,j}.inf);
                J2(i,j) = interval(Fmin(i,j),Fmax(i,j));
            end
        end
        J = J2;
    otherwise
        disp('Wrong Argument');
end

switch (h_str)
    case 'C1'
        h = center(interval_X);
    case 'C2'
        h = h_C2(X0,J,1);
    otherwise
        disp('Wrong Argument');
end

Z = func(h,center(interval_W));
Z = conZonotope(Z);

% Predicted Set
X_predicted = Z + theorem_01(J,X0+(-h));

end

function S = theorem_01(J,X)

% Constraint Reduction (X0 -> Xbar)
if isempty(X.A)
    Xbar = X;
else
    Xbar = reduce(X,'girard',1,0);
end

pbar = Xbar.Z(:,1);     % Center of Reduced Zonotope
Mbar = Xbar.Z(:,2:end); % Generators of Reduced Zonotope

% P Matrix
ngbar = size(Mbar,2);
n = size(J,2);

m = (J-center(J))*pbar;

P = zeros(1,n);
for i = 1:n
    term = 0;
    for j = 1:ngbar
        for k = 1:length(m)
            term = term+((2*rad(J(i,k)*abs(Mbar(k,j)))));
        end
    end
    P(i) = 1/2*(2*rad(m(i))) + 1/2*term;
end

P = diag(P);

% S (Output of Theorem 01) -> S = center(J)*X + P B_inf^n
S1  = center(J)*X;

Sc = S1.Z(:,1);
SG = [S1.Z(:,2:end) P];

if isempty(S1.A)
    S = conZonotope(Sc,SG);
else
    SA = [S1.A zeros(size(S1.A,1),n)];
    Sb = S1.b;
    S = conZonotope(Sc,SG,SA,Sb);
end

end