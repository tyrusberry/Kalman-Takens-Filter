function [stateEstimate,varEstimate,Q,R] = AUKF(state,obs,f,h,Q_param,stationarity,Q,R)

% state         an Nx1 vector, best estimate of the initial state
% obs           an MxT matrx which contains T observations in R^M
% f             a function handle which computes the dynamics as a 
%                function of the state
% h             a function handle which computes the observation as a 
%                function of the state
% Q_param       0 - none, 1 - full, 2 - multiple of identity, 3 - diagonal
% stationarity  widow length for adaptively averaging Q and R
% Q             the NxN covariance matrix for the dynamical noise
% R             the MxM covariance matrix for the observation noise

    N = size(state,1);
    M = size(obs,1);
    T = size(obs,2);

    if (nargin<8) R = eye(M); end
    if (nargin<7) Q = eye(N); end
    if (nargin<6) stationarity = 2000; end
    if (nargin<5) Q_param = 1; end

    covmat = eye(N);
    
    stateEstimate = zeros(N,T);
    varEstimate = zeros(N,T);
    prev_innov = zeros(M,1);
    prev_FPF = eye(N);

    stateEstimate(:,1) = state;
    
    for t = 1:T-1

        [state,covmat,Qest,Rest,prev_innov,prev_FPF] = AUKF_STEP(state,covmat,obs(:,t+1),f,h,Q,R,Q_param,prev_innov,prev_FPF);         
        
        stateEstimate(:,t+1) = state;
        varEstimate(:,t+1) = diag(covmat);
        
        if (Q_param > 0)               
            if (t>20)
                Q = (Q*(2*stationarity-1) + Qest)/(2*stationarity);
                R = (R*(stationarity-1) + Rest)/(stationarity);
            end   
        end
        
    end

    %%% fix up the final estimates of Q and R
    [U,S,~] = svd(Q);
    Q = U*S*U';
    [U,S,~] = svd(R);
    R = U*S*U';
    
