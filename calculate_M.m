function [MM] = calculate_M(locMap,triMap,phi)
%Efficient approach (sparse matrix)
M = length(locMap);
T = length(triMap);

i = zeros(9*T, 1);
j = zeros(9*T, 1);
mm = zeros(9*T, 1);

phi_mdpt = phi(2);

count = 0;
for R = 1:1:T
    %calculate area of triangle
    rix = locMap( triMap(R,1), 1 );
    riy = locMap( triMap(R,1), 2 );
    rjx = locMap( triMap(R,2), 1 );
    rjy = locMap( triMap(R,2), 2 );
    rkx = locMap( triMap(R,3), 1 );
    rky = locMap( triMap(R,3), 2 );

    areaR = abs(det([1, 1, 1; rix, rjx, rkx; riy, rjy, rky]))/2;

    %calculate local Mass matrix
    if phi_mdpt == 0
        A = [1 0 0;
             0 1 0;
             0 0 1];
    else
        A = [2 1 1;
             1 2 1;
             1 1 2] * phi_mdpt^2;
    end
    MMlocal = areaR/3 * A;
    
    %sparse elements
    for ti = 1:3
        for tj = 1:3
            count = count+1;
            i(count) = triMap(R,ti);
            j(count) = triMap(R,tj);
            mm(count) = MMlocal(ti,tj);
        end;
    end;
end;

MM = sparse(i,j,mm,M,M);

end


%%Inefficient but more natural approach
%M = length(locMap);
%MM = zeros(M,M);
%
%for R = 1:1:length(triMap)
%    %calculate area of triangle
%    rix = locMap( triMap(R,1), 1 );
%    riy = locMap( triMap(R,1), 2 );
%    rjx = locMap( triMap(R,2), 1 );
%    rjy = locMap( triMap(R,2), 2 );
%    rkx = locMap( triMap(R,3), 1 );
%    rky = locMap( triMap(R,3), 2 );
%
%    areaR = abs(det([1, 1, 1; rix, rjx, rkx; riy, rjy, rky]))/2;
%
%    %calculate local Mass matrix
%    MMlocal = areaR * [1/6 1/12 1/12; 1/12 1/6 1/12; 1/12 1/12 1/6];
%
%    %sum into Mass matrix
%    for i = 1:3
%        for j = 1:3
%            MM(triMap(R,i),triMap(R,j)) = MM(triMap(R,i),triMap(R,j)) + MMlocal(i,j);
%        end;
%    end;
%
%end;
