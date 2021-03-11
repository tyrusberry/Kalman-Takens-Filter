clear;clc;close all;
addpath('AdaptiveUKF');
addpath('KalmanTakens');
addpath('Lorenz96');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 40;                         %%% number of state variables
T = 5000;                       %%% total number of time steps
dt = 0.05;                      %%% time between observations
dynnoiseVariance = 0;           %%% system noise (stochastic forcing)
obsnoiseVariance = 16;          %%% observation noise

Q = dynnoiseVariance*eye(N);    %%% system noise covariance matrix
R = obsnoiseVariance*eye(N);    %%% obs noise covariance matrix

[truth,obs,p] = GenerateL96(T,N,dt,Q,R);


%%%%%%%%%%%%%%%%%% Run the filter with the true dynamics %%%%%%%%%%%%%%%%%%

f=@(x) L96Dynamics(x,dt,p);     %%% true dynamics
h=@(x) L96Obs(x,N);             %%% true obs function
state = obs(:,1);               %%% initial state

tic;
stateEstimateTrueModel = AUKF(state,obs,f,h);    %%% Run the filter
toc;


%%%%%%%%%%%%%%% Run the Kalman-Takens filter with no model %%%%%%%%%%%%%%%%

tic;
for i = 1:N                     %%% run a filter for each observation

    %%% It sometimes helps to include spatially neighboring states
    timeSeries =[obs(i,:); obs(rema(i-1,N),:); obs(rema(i+1,N),:)];
    %timeSeries = obs(i,:);     %%% this would only use the current state
    
    delays = 4;         %%% number of delays in the Takens embedding
    kNN=40;             %%% number of nearest neighbors to use in forecast
    
    %%% Run the Kalman-Takens Filter which uses the same UKF code as above
    %%% but with no model
    stateEstimatei = KalmanTakensFilter(timeSeries,delays,kNN);
    
    % the first component of the K-T output is the i-th state estimate
    stateEstimateNoModel(i,:)=stateEstimatei(1,:);      
    
end
toc;


%%%%%%%%%%%%%%%%%%%% Compute errors and plot results %%%%%%%%%%%%%%%%%%%%%%

filterTransient=2000;   %%% use last 3000 points to compute errors
obsRMSE = sqrt(mean(mean((obs(:,filterTransient:end)-truth(:,filterTransient:end)).^2)))
trueModelRMSE = sqrt(mean(mean((stateEstimateTrueModel(:,filterTransient:end)-truth(:,filterTransient:end)).^2)))
noModelRMSE = sqrt(mean(mean((stateEstimateNoModel(:,filterTransient:end)-truth(:,filterTransient:end)).^2)))

plotInds = T-1000:T;    %%% plot the last 1000 points

figure(1);
plot(truth(1,plotInds),'k','linewidth',2);hold on;
plot(obs(1,plotInds),'go');
plot(stateEstimateTrueModel(1,plotInds),'r');
plot(stateEstimateNoModel(1,plotInds),'b');
legend('Truth','Observations','True Model Filter','No Model Filter');
xlim([0 1000]);

figure(2);
subplot(4,1,1);imagesc(truth(:,plotInds),[-10 15]);title('Truth')
subplot(4,1,2);imagesc(obs(:,plotInds),[-10 15]);title('Observations')
subplot(4,1,3);imagesc(stateEstimateTrueModel(:,plotInds),[-10 15]);title('True Model Filter Estimate')
subplot(4,1,4);imagesc(stateEstimateNoModel(:,plotInds),[-10 15]);title('No Model Filter Estimate')



    
