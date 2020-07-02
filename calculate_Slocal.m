function [Slocal] = calculate_Slocal(R,locMap,triMap)

Slocal = zeros(3,3);

%area of triangle
rix = locMap( triMap(R,1), 1 );
riy = locMap( triMap(R,1), 2 );
rjx = locMap( triMap(R,2), 1 );
rjy = locMap( triMap(R,2), 2 );
rkx = locMap( triMap(R,3), 1 );
rky = locMap( triMap(R,3), 2 );

areaR = abs(det([1, 1, 1; rix, rjx, rkx; riy, rjy, rky]))/2;

%transform Af
Af = [rix-rkx, riy-rky; rjx-rkx, rjy-rky];

%gradient function
phi = [[1,0]', [0,1]', [-1,-1]'];

%local Stiffness matrix
for i = 1:3
    for j = 1:3
        Slocal(i,j) = areaR * ( (Af\phi(:,i))' * (Af\phi(:,j)) );
    end;
end;


