function [dJ,dcon,dconeq] = Backward_Implicit(theta,a,b,rho_obs,n,v,f,phi,lambda0,c,dt,N,M,w)

kappa = theta(1);
D = theta(2);
lambda1 = theta(3);
mu = theta(end);

%Step 0a: Calculate Mass and Stiffness matrices
MM = calculate_M(v,f,phi);
S = calculate_S(v,f);

%Step 0b: Initialize Intermediary Variables
onesm = ones(1,M);
z = zeros(N,M);
X = full((MM + dt * D * S)');

%Step 1: Initialize pN = gradient g(xN)
pa = zeros(N,M);
pb = zeros(N,M);

%pa(N,:) = 0;
for i = 1:1:length(n)
    rho_max_aged = (mu*c) .* exp(lambda0*dt*(double(n(i))-1));
    pb(n(i),:) = -2 * (rho_max_aged ./ (onesm + exp(b(n(i),:))) - rho_obs(i,:)) .* (rho_max_aged .* exp(b(n(i),:)) ./  (onesm + exp(b(n(i),:))).^2) ;
end

%Step 2: Recursively Calculate p
for i = N:-1:2 %tmax:dt:tmin
    z(i,:) = ( linsolve(X, pa(i,:)') )'; 
    dF = calculate_dF(a(i-1,:),v,f, phi);
    [panew, pbnew] = backward_implicit(z(i,:), pb(i-1,:), pb(i,:), dt, kappa, lambda1, MM, dF);
    pa(i-1,:) = panew;
    pb(i-1,:) = pbnew;
end
z(1,:) = ( linsolve(X, pa(1,:)') )'; 

%Step 3: Update dF
%kappa, D
dJ(1) = 0;
dJ(2) = 0;
for i = 1:1:N
    F = calculate_F(a(i,:)',v, f, phi); 
    dJ(1) = dJ(1) + z(i,:) * dt * F;
    dJ(2) = dJ(2) +  -z(i,:) * (dt * S) * ((MM + dt * D * S)\(MM*a(i,:)' + dt * kappa * F));
end

%lambda1
dJ(3) = dt * sum(sum(pb .* a));

%a0
dF = calculate_dF(a(1,:),v,f, phi);
da0J = z(1,:)*(MM + dt*kappa*dF) + pb(1,:)*lambda1*dt;
sigma = 1;
dcxa0 = (v(:,1)-theta(4)*onesm')./sigma^2 .* a(1,:)';
dcya0 = (v(:,2)-theta(5)*onesm')./sigma^2 .* a(1,:)';
dJ(4) = da0J*dcxa0;
dJ(5) = da0J*dcya0;

dcon =   [0 0 0 -w*dcxa0 -w*dcya0 0 0];
dconeq = [0 0 0 0 0 0 0];

%b0
dJ(6) = sum(pb(1,:));

%mu
for i = 1:length(n)
    numerator(i,:) = mu*onesm ;
    denominator(i,:) = c .* exp(lambda0*dt*(double(n(i))-1)) ./ ( onesm +exp(b(int16(n(i)),:)) );
end
dJ(7) = 2* sum(sum( numerator.*(denominator.^2) - rho_obs.*denominator ));

end
