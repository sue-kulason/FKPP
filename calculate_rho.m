function [rho] = calculate_rho(theta,b,n,lambda0,dt,M)
for i = 1:length(n)
    rho(i,:) = (theta(end-M+1:end)) .* exp(lambda0*dt*(double(n(i))-1)) ./ (ones(1,M)+exp(b(n(i),:)));
end