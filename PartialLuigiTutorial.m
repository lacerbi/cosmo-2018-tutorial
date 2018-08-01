%% CoSMo 2017, Day 5 -- Model fitting & optimization tutorial with Luigi Acerbi

% Prepare data, from Konrad's tutorial

% Load data
load M1_Stevenson_Binned

% Remove all times where speeds are very slow
isGood=find(handVel(1,:).^2+handVel(2,:).^2>.015);
handVel=handVel(1:2,isGood);
handPos=handPos(1:2,isGood);
spikes=spikes(:,isGood);
time=time(isGood);
angle=atan2(handVel(1,:),handVel(2,:));

% Plot Raw Data
nNeuron = 193  %193
clf
hold on
plot(angle,spikes(nNeuron,:)+0.2*randn(size(spikes(nNeuron,:))),'r.')

%% Part I: Model fitting
% Fix random seed for reproducibility
rng(1);

% Define exponential-cosine tuning curve

% Define optimization starting point and bounds

% Randomize initial starting point inside plausible box

%% Fit with (bounded) fminsearch first


%% Fit with fmincon (may get stuck!)


%% Fit with BADS (download from: https://github.com/lacerbi/bads)



