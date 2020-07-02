function [con,coneq] = myconstraint(theta,v,w)

sigma = 1;
A = sqrt(2*pi)*sigma;

a_0 = theta(4:5);
[M,~] = size(v);
a = initialize_a(1,M,a_0(1),a_0(2),v);
con = .95*A - w*a';
coneq = 0;

[theta(1:6) mean(theta(7:end))]
con