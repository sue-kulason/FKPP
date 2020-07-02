function [J,dJ,dcon,dconeq] = fkppwithgrad(theta,rho_obs,t,v,f,phi,lambda0,c,dt,w)

tmax = max(t);
tmin = 50;
n = int16( (t-tmin)/dt +1 );
N = size(tmin:dt:tmax,2);
M = size(v,1);

[a,b] = Forward_Implicit(theta,v,f,phi,dt,N,M);
rho = calculate_rho(theta,b,n,lambda0,c,dt);
J = sum(sum( (rho - rho_obs).^2 ));
[dJ,dcon,dconeq] = Backward_Implicit(theta,a,b,rho_obs,n,v,f,phi,lambda0,c,dt,N,M,w);

