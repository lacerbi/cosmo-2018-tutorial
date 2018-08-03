%% CoSMo 2018, Day 5 -- Model fitting & optimization tutorial with Luigi Acerbi

%% Add utilities' folder to MATLAB path

baseFolder = fileparts(which('PartialLuigiTutorial.m'));
addpath([baseFolder,filesep(),'utils']);

%% Prepare Stevenson et al.'s (2011) data, from Konrad's tutorial

% Load data into a struct
data = load('M1_Stevenson_Binned.mat');

% Remove all times where speeds are very slow
isGood=find(data.handVel(1,:).^2+data.handVel(2,:).^2>.015);
data.handVel=data.handVel(1:2,isGood);
data.handPos=data.handPos(1:2,isGood);
spikes_vec=data.spikes(:,isGood);
data.time=data.time(isGood);
angle_vec=atan2(data.handVel(1,:),data.handVel(2,:));

% Plot Raw Data
nNeuron = 193  %193
clf
hold on
plot(angle_vec,spikes_vec(nNeuron,:)+0.2*randn(size(spikes_vec(nNeuron,:))),'r.')
xlabel('Angle (rad)');  xlim([-pi,pi]);
ylabel('Spikes (Hz)');
set(gca,'TickDir','out','Box','off');
set(gcf,'Color','w');

%% Model fitting!

% Fix random seed for reproducibility
rng(1);

% First, write the negative log likelihood (nLL) function in a separate file

% Define exponential-cosine tuning curve
tuningFun = @(params,angle) exp(params(1)+params(2).*cos(angle-params(3)));

% Alternative tuning curve, rectified cosine
% tuningFun = @(params,angle) max(params(1)+params(2).*cos(angle-params(3)),1e-4);   % Why did I put 1e-4 and not 0?

% Define handle to negative log likelihood function
nLLfun = @(params) Tuning_nLL(params, spikes_vec(nNeuron,:), angle_vec, tuningFun);

% Define optimization starting point and bounds
LB = [-10,0.001,-pi];     % Lower bound
UB = [50,20,pi];  % Upper bound
PLB = [0,0.1,-pi];     % Plausible lower bound
PUB = [5,2,pi];   % Plausible upper bound

% Randomize initial starting point inside plausible box
x0 = rand(size(PLB)).*(PUB - PLB) + PLB;


%% Fit with (bounded) fminsearch first
options.Display = 'iter';
options.OutputFcn = @(x, optimValues, state) plotProfile(nLLfun,x,PLB,PUB,[],{'b_0','b_1','\theta_0'});

%[bestParams,fvalCosExp] = fminsearchbnd(nLLfun,x0,LB,UB,options);
% plot(-pi:pi/80:pi, exp(bestParams(1)+bestParams(2)*cos((-pi:pi/80:pi)-bestParams(3))))
%bestParams
%fvalCosExp


%% Fit with fmincon (may get stuck!)
% [bestParams,fvalCosExp] = fmincon(nLLfun,x0,[],[],[],[],LB,UB,[],options);
% plot(-pi:pi/80:pi, exp(bestParams(1)+bestParams(2)*cos((-pi:pi/80:pi)-bestParams(3))))
% bestParams
% fvalCosExp

%% Fit with BADS (download from: https://github.com/lacerbi/bads)

options.Display = 'iter';
options.UncertaintyHandling = false;    % Log likelihood is exact
options.PeriodicVars = 3;               % Preferred direction is periodic

[bestParams,fvalCosExp,exitflag,output] = bads(nLLfun,x0,LB,UB,PLB,PUB,[],options);
% plot(-pi:pi/80:pi, exp(bestParams(1)+bestParams(2)*cos((-pi:pi/80:pi)-bestParams(3))))
bestParams
fvalCosExp
output

%% Fit with CMA-ES (very powerful, but uses *a lot* of function evaluations)
[bestParams,fvalCosExp,~,output] = cmaes_wrapper('Tuning_nLL',x0,LB,UB,PLB,PUB,[],spikes_vec(nNeuron,:),angle_vec,tuningFun);
bestParams
fvalCosExp
output

%% VBMC special preview! (to appear at: https://github.com/lacerbi/vbmc)

rng(1);
%x0 = rand(size(PLB)).*(PUB - PLB) + PLB;
x0 = bestParams;    % Better to initialize from known optimum
vp = vbmc(@(x) -nLLfun(x),x0,LB,UB,PLB,PUB,struct('Display','iter','Plot','on'));
pause


%% MCMC via slice sampling

N = 2e3;
smpl_options = [];

[samples,fvals,exitflag,output] = slicesamplebnd(@(x) -nLLfun(x),x0,N,[],LB,UB,smpl_options);
cornerplot(samples,{'b_0','b_1','\theta_0'},[],[PLB;PUB]);

%% Part II: Fitting Drifting Tuning Curves

rng(1);
npivots = 10;

% Define optimization starting point and bounds
LB = [-20,0,0,-2*pi*ones(1,npivots-1)];     % Lower bound
UB = [20,10,2*pi,2*pi*ones(1,npivots-1)];  % Upper bound
PLB = [-3,0,0,-pi*ones(1,npivots-1)];     % Plausible lower bound
PUB = [3,3,2*pi,pi*ones(1,npivots-1)];   % Plausible upper bound

x0 = rand(size(PLB)).*(PUB - PLB) + PLB;

nLLfun = @(params) DriftingTuning_nLL(params, spikes_vec(nNeuron,:), angle_vec, time, tuningFun);

% [bestParams,fvalCosExp] = fminsearchbnd(nLLfun,x0,LB,UB,options);

% Fit with BADS (https://github.com/lacerbi/bads)

options.Display = 'iter';
options.UncertaintyHandling = false;    % Log likelihood is exact
options.PeriodicVars = 3;       % Preferred directions are periodic
options.OutputFcn = [];
% @(x, optimValues, state) plotProfile(nLLfun,x,PLB,PUB);

[bestParams,fvalCosExp,exitflag,output] = bads(nLLfun,x0,LB,UB,PLB,PUB,[],options);
%plot(-pi:pi/80:pi, exp(bestParams(1)+bestParams(2)*cos((-pi:pi/80:pi)-bestParams(3))))
fvalCosExp

%% MCMC via slice sampling

N = 8e3;
smpl_options = [];

[samples,fvals,exitflag,output] = slicesamplebnd(@(x) -nLLfun(x),x0,N,[],LB,UB,smpl_options);
cornerplot(samples,[],[],[PLB;PUB]);

%% Cart-pole inverse dynamics

% Add Konrad's folder
baseFolder = fileparts(which('PartialLuigiTutorial.m'));
addpath([baseFolder,filesep(),'Konrad']);

m1=1;
m2=1;
g=9.81;
l=0.3;
deltaT=0.0001;
duration=20;

X=[0 pi-pi/16]; %pi=up, 0/2pi=down
dotX=[0 0];
x0=[1 2 0.01];
LB = [-40,-40,-20];
UB = [40,40,20];
PLB = [-10,-10,-10];
PUB = [10,10,10];
options.Display = 'iter';
lossFun = @(x) brutalScore(x,X,dotX,m1,m2,g,l);

options.OutputFcn = @(x, optimValues, state) plotProfile(lossFun,x,LB,UB,64);
%[bestParas, score]=bads(lossFun,x0,LB,UB,[],[],options);

% [bestParas, score]=fminsearchbnd(lossFun,x0,LB,UB,options);

[bestParams,fval,exitflag,output] = bads(lossFun,x0,LB,UB,PLB,PUB,[],options);
bestParams
fval

[bestParams,fval,~,output] = cmaes_wrapper('brutalScore',x0,LB,UB,PLB,PUB,[],X,dotX,m1,m2,g,l);
bestParams
fval
