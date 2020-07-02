function [F] = calculate_F(a,locMap, triMap, phi)

M = length(locMap);
T = length(triMap);

i = zeros(3*T, 1);
j = zeros(3*T, 1);
f = zeros(3*T, 1);

count = 0;
for R = 1:1:T

    %calculate local Stiffness matrix
    Flocal = calculate_Flocal(a, R, locMap, triMap, phi);

    %sum into Stifness matrix
    for ti = 1:3
        count = count+1;
        i(count) = triMap(R,ti);
        j(count) = 1;
        f(count) = Flocal(ti);
    end;

end;

F = sparse(i,j,f,M,1);

end