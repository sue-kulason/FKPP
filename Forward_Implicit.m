function [a,b] = Forward_Implicit(theta,v,f,phi,dt,N,M)

MM = calculate_M(v,f,phi);
S = calculate_S(v,f);

kappa = theta(1);
D = theta(2);
lambda1 = theta(3);
a_0 = theta(4:5);
b_0 = theta(6);

%Initialize Variables
a = initialize_a(N,M,a_0(1),a_0(2),v);
b = initialize_b(N,M,b_0);

%Calculate x
for i = 1:1:N-1
  F = calculate_F(a(i,:)',v, f, phi);
  [anew, bnew] = forward_implicit(a(i,:)', b(i,:)', dt, kappa, D, lambda1, MM, S, F);
  a(i+1,:) = anew';
  b(i+1,:) = bnew';
end