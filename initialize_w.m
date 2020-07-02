function [w] = initialize_w(v,f,phi)

[M,~] = size(v);
[T,~] = size(f);
w = zeros(1,M);

for R = 1:1:T
    %calculate area of triangle
    rix = v( f(R,1), 1 );
    riy = v( f(R,1), 2 );
    rjx = v( f(R,2), 1 );
    rjy = v( f(R,2), 2 );
    rkx = v( f(R,3), 1 );
    rky = v( f(R,3), 2 );

    areaR = abs(det([1, 1, 1; rix, rjx, rkx; riy, rjy, rky]))/2;
    
    w(f(R,1)) = w(f(R,1)) + phi(end)*areaR;
    w(f(R,2)) = w(f(R,2)) + phi(end)*areaR;
    w(f(R,3)) = w(f(R,3)) + phi(end)*areaR;
end