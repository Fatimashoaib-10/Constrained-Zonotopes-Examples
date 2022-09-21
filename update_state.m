function Xhat = update_state(Xbar,yk,C,Du,uk,Dv,V)

YK = (yk-Du*uk)+(-Dv*V);

Xhat = cz_intersection(Xbar,YK,C);

end

function res = cz_intersection(cz1,cz2,R)

% Add trivial constraints if not already present

if isempty(cz1.A)
    cz1A = zeros(1,size(cz1.Z,2)-1);
    cz1b = 0;
    cz1 = conZonotope(cz1.Z,cz1A,cz1b);
end

if isempty(cz2.A)
    cz2A = zeros(1,size(cz2.Z,2)-1);
    cz2b = 0;
    cz2 = conZonotope(cz2.Z,cz2A,cz2b);
end

% Intersection

Z = [cz1.Z,zeros(size(cz2.Z)-[0,1])];
A = blkdiag(cz1.A,cz2.A);
A = [A; R*cz1.Z(:,2:end), -cz2.Z(:,2:end)];
b = [cz1.b; cz2.b; cz2.Z(:,1)-R*cz1.Z(:,1)];

res = conZonotope(Z,A,b);

end