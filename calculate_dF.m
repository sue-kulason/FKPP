function [dF] = calculate_dF(a,locMap, triMap, phi)

M = length(locMap);
T = length(triMap);

i = zeros(3*T, 1);
j = zeros(3*T, 1);
df = zeros(3*T, 1);

count = 0;
for R = 1:1:T

    %calculate local Stiffness matrix
    dFlocal = calculate_dFlocal(a, R, locMap, triMap, phi);

    %sum into Stifness matrix
    for ti = 1:3
        for tj = 1:3
            count = count+1;
            i(count) = triMap(R,ti);
            j(count) = triMap(R,tj);
            df(count) = dFlocal(ti,tj);
        end
    end;

end;

dF = sparse(i,j,df,M,M);
