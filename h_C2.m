function h = h_C2(X,J,nc)

eta = sdpvar(size(X.Z(:,2:end),2),1);

n = size(J,2);

%% O
O = zeros(1,n);
for i = 1:n
    tem = 0;
    for j = 1:n
        tem = 2*rad(J(i,j));
    end
    O(i) = tem;
end

O = diag(O);

%% p_bar

term = 0;
Xcurrent  = X;
Xpre      = cell(1,nc);
Xrescaled = cell(1,nc);
Xreduced  = cell(1,nc);
etam      = cell(1,nc);
AG        = cell(1,nc);

for i = 1:nc

Xpre{i}      = precond(Xcurrent);             % Preconditioning
Xrescaled{i} = rescale(Xcurrent,'iter');      % Rescaling
Xreduced{i}  = reduce(Xcurrent,'girard',0,0); % Constraint Reduction

% etam{i} = Xpre{i}.Z(:,2:end)\(Xrescaled{i}.Z(:,1) - Xcurrent.Z(:,1));
etam{i} = eta_m(Xpre{i});

% AG{i} = Xrescaled{i}.b\(Xreduced{i}.Z(:,1) - Xrescaled{i}.Z(:,1));
AG{i} = A_G(Xrescaled{i},size(Xrescaled{i}.A,2));

term = term + Xpre{i}.Z(:,2:end)*etam{i} + AG{i}*Xrescaled{i}.b;

Xcurrent = Xreduced;
end

Gx = X.Z(:,2:end);
p_bar = -Gx*eta + term;

%%
constr  = [X.A*eta == X.b, norm(eta,inf)<=1 ];%-1<=eta<=1];

options = sdpsettings('verbose',0);
obj     = norm(O*p_bar,1);
optimize(constr, obj, options); 

h = X.Z(:,1) + Gx*value(eta);

end

%% Functions %% (Modified - Original from CORA 2020)

function Xpre = precond(X)
% preconditioning of the constraint matrix to improve rescaling    

    % bring constraint matrices to Reduced Echelon Form
    [A_,b_,indPer] = gauss_jordan(X.A,X.b);
    
    % adapt the generator matrix to the new order of factors
    c = X.Z(:,1);
    G = X.Z(:,2:end);
    G = G(:,indPer);
   
    Xpre = conZonotope(c,G,A_,b_);
end

function [A_,b_,indPer] = gauss_jordan(A,b)
% Transform the matrix [A b] to Reduced Echelon Form using Gauss-Jordan
% elimination with full pivoting. The row elements with the largest 
% absolute values relative to the inifinity norm of their row are chosen as 
% the pivot elements.

    A_ = [A,b];
    [m,n] = size(A_);
    
    tol = max(m,n)*eps(class(A_))*norm(A_,'inf');

    % loop over the entire matrix
    indPer = 1:n-1;

    for i = 1:m
        
       % unreduced submatrix
       Atemp = A_(i:m,i:n-1);
        
       % divide each row by it's inifinity norm
       normInf = sum(abs(Atemp),2);
       Atemp = diag(1./normInf) * Atemp;
        
       % find value and index of largest element in the remainder of column j
       if size(Atemp,1) > 1
           [valTemp,indTemp] = max(abs(Atemp));
           [p,indC] = max(valTemp);
           indR = indTemp(indC) + i - 1;
           indC = indC + i - 1;
       else
          [p,indC] = max(abs(Atemp));
          indC = indC + i - 1;
          indR = i;
       end
       
       if (p > tol)
          
          % bring the row with the pivot element up
          A_([i indR],:) = A_([indR i],:);
          
          % bring the column with the pivot element to the front
          A_(:,[i indC]) = A_(:,[indC i]);
          temp = indPer(i);
          indPer(i) = indPer(indC);
          indPer(indC) = temp;
          
          % divide the pivot row by the pivot element
          Ai = A_(i,:)/A_(i,i);  
          
          % subtract multiples of the pivot row from all the other rows
          A_(:,i:n) = A_(:,i:n) - A_(:,i)*Ai(:,i:n);
          A_(i,i:n) = Ai(i:n);
       end
    end
    
    b_ = A_(:,end);
    A_ = A_(:,1:n-1);
end

function ksi_m = eta_m(obj)

if isempty(obj.ksi)
    obj = precond(obj);
    domKsi = eta_iterative(obj);
    ksi_l = infimum(domKsi);
    ksi_u = supremum(domKsi);
else
    ksi_l = min(obj.ksi,[],2);
    ksi_u = max(obj.ksi,[],2);
end

ksi_m = (ksi_u + ksi_l)/2;

end

function Dksi = eta_iterative(obj)

iter = 1;

[m,n] = size(obj.A);
E = interval(-ones(n,1), ones(n,1));
R = interval(-Inf(n,1), Inf(n,1));

A = obj.A;
b = obj.b;
iA = A.^-1;

for k = 1:iter
    for i = 1:m
        for j = 1:n
            if iA(i,j) ~= Inf
                temp = E;
                temp(j) = 0;
                dummy = iA(i,j) .* ( b(i) - A(i,:)*temp );
                R(j) = R(j) & dummy;
                E(j) = E(j) & R(j);
            end
        end
    end
end
Dksi = E;

end

function AG = A_G(obj,ind)

G = obj.Z(:,2:end);
A = obj.A;

[m,n] = size(A);
ind1 = find(A(:,ind) ~= 0);

if ~isempty(ind1)
    ind1 = ind1(1);
    a = A(ind1, ind);
    E = zeros(n,m);
    E(ind, ind1) = 1;
    
    AG = G * E ./a;
else
    AG = [];
end

end