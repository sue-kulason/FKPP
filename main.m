% User Input
subject = 'example/test1.mat';
initialguess = 'example/initialguess3.mat';
resultDir = 'example/';

dt = 0.1;
phi = [1, 1/2, 1/3];

% Load in Data
load(subject);
temp = strsplit(subject,{'.','/'});
subj = temp{end-1};

load(initialguess);
temp = strsplit(initialguess,{'.','_','/'});
idx = temp{end-1};

w = initialize_w(v,f,phi);

% Perform Gradient Descent
objfunc = @(theta) fkppwithgrad(theta,rho,t,v,f,phi,lambda0,dt,w);
confunc = @(theta) myconstraint(theta,v,w);

options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                       'MaxFunctionEvaluations',3000,'MaxIterations',1000,... 
                       'Display','iter-detailed');
M = size(v,1);
lb = [.1, .1, .01, min(v), -10, 1*ones(1,M)];
ub = [10, 10,   1, max(v),  -1, 6*ones(1,M)];

tic;
[theta,J,flag,output] = fmincon(objfunc,theta0,[],[],[],[],lb,ub,confunc,options);
time = toc;

% Save Results
tmax = max(t);
tmin = 50;
n = int16( (t-tmin)/dt +1 );
N = size(tmin:dt:tmax,2);
M = size(v,1);
[a,b] = Forward_Implicit(theta,v,f,phi,dt,N,M);
rho_model = calculate_rho(theta,b,1:N,lambda0,dt,M);

save([resultDir '/' subj '_' idx '_output.mat'])
