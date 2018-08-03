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

% Define optimization starting point and bounds

% Randomize initial starting point inside plausible box

%% Fit with (bounded) fminsearch first


%% Fit with fmincon (may get stuck!)

%% Fit with CMA-ES (very powerful, but uses *a lot* of function evaluations)



%% Fit with BADS (download from: https://github.com/lacerbi/bads)

%% Another problem: cart-pole inverse dynamics

% % Add Konrad's folder
% baseFolder = fileparts(which('PartialLuigiTutorial.m'));
% addpath([baseFolder,filesep(),'Konrad']);
% 
% m1=1;
% m2=1;
% g=9.81;
% l=0.3;
% deltaT=0.0001;
% duration=20;
% 
% X=[0 pi-pi/16]; %pi=up, 0/2pi=down
% dotX=[0 0];
% x0=[1 2 0.01];
% LB = [-40,-40,-20];
% UB = [40,40,20];
% PLB = [-10,-10,-10];
% PUB = [10,10,10];
% options.Display = 'iter';
% lossFun = @(x) brutalScore(x,X,dotX,m1,m2,g,l);
% 
% options.OutputFcn = @(x, optimValues, state) plotProfile(lossFun,x,LB,UB,64);
% 
% % Try and fit it with different optimizers!
