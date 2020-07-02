function [Flocal] = calculate_Flocal(a,R,locMap, triMap,phi)

%area of triangle
rix = locMap( triMap(R,1), 1 );
riy = locMap( triMap(R,1), 2 );
rjx = locMap( triMap(R,2), 1 );
rjy = locMap( triMap(R,2), 2 );
rkx = locMap( triMap(R,3), 1 );
rky = locMap( triMap(R,3), 2 );

areaR = abs(det([1, 1, 1; rix, rjx, rkx; riy, rjy, rky]))/2;

%note that A is symmetric
A = [ [a(triMap(R,1)), a(triMap(R,2)), a(triMap(R,3))]', ...
      [a(triMap(R,2)), a(triMap(R,3)), a(triMap(R,1))]', ...
      [a(triMap(R,3)), a(triMap(R,1)), a(triMap(R,2))]' ]; 

B = phi * [3 0 0;
           8 4 4;
           9 9 9];

C = [phi * [3 0 0;
            4 2 2;
            3 3 3];
     phi * [0 0 0;
            2 2 0;
            3 3 3];
     phi * [0 0 0;
            2 0 2;
            3 3 3];
     ];

%local F matrix
Flocal = (areaR/60) * [B*A(:,1)-A(:,1)'*C*A(:,1); ...
                       B*A(:,2)-A(:,2)'*C*A(:,2); ...
                       B*A(:,3)-A(:,3)'*C*A(:,3)];
end
