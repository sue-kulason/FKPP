function [S] = calculate_S(locMap,triMap)

M = length(locMap);
T = length(triMap);

i = zeros(9*T, 1);
j = zeros(9*T, 1);
s = zeros(9*T, 1);

count = 0;
for R = 1:1:T

    %calculate local Stiffness matrix
    Slocal = calculate_Slocal(R, locMap, triMap);

    %sum into Stifness matrix
    for ti = 1:3
        for tj = 1:3
            count = count+1;
            i(count) = triMap(R,ti);
            j(count) = triMap(R,tj);
            s(count) = Slocal(ti,tj);
        end;
    end;

end;

S = sparse(i,j,s,M,M);

