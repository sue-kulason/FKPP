function a = initialize_a(N,M,x,y,locMap)

ctr = [x,y];
sigma = 1; 
a = zeros(N,M);

for i = 1:M
    dist = locMap(i,:) -ctr;
    a(1,i) = exp(-(dist*dist')/(2*sigma^2))/sqrt(2*pi*sigma^2);
end

end