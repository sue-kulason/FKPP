function [rho] = calculate_rho(theta,b,n,lambda0,c,dt)
for i = 1:length(n)
    rho(i,:) = (theta(end).*c) .* exp(lambda0*dt*(double(n(i))-1)) ./ (ones(size(c))+exp(b(n(i),:)));
end