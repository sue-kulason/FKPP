function [a,b] = forward_implicit(a,b,dt,kappa,D,lambda1,MM,S,F)

b = b + dt*lambda1*a;
a = (MM + dt * D * S) \ (MM * a) + (MM + dt * D * S) \ ( dt * kappa * F);
a(a<0) = 0;
a(a>1) = 1;

end
