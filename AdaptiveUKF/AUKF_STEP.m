function [xa,Pa,Qest,Rest,innov,FPF] = AUKF_STEP(xa,Pa,obs,f,h,Q,R,Q_param,prev_innov,prev_FPF) 
%%%Propagate from one observation time step to the next%%%

% xa        an Nx1 vector containing the state at itime
% Pa        an NxN matrix containing the covariance matrix at itime
% obs      an Mx1 vector which contains the observation at ftime
% f         a function handle which computes the vector field as a 
%            function of the state and time
% h         a function handle which computes the observation as a 
%            function of the state and time
% Q         the NxN covariance matrix for the dynamical noise
% R         the MxM covariance matrix for the observation noise
% Q_param   0 - none, 1 - full, 2 - multiple of identity, 3 - diagonal


warning('off','MATLAB:rankDeficientMatrix');

    N = length(xa);
    M = length(obs);

    %%% Fix up the covariance matrices to insure symmetric and pos. def. %%%
    
    [U,S,~] = svd(Q);
    Q = U*S*U';
    [U,S,~] = svd(R);
    R = U*S*U';
    
    %%% Generate an unscented ensemble, X, and associated weights %%%
    
    [X,weights] = unscentedEnsemble(xa,Pa,sqrt(N-1));
    weightMatrix = diag(weights(2:end));

    %%% Apply the nonlinear dynamics, f, and observation, h, to the ensemble members
    
    FX = f(X);
    Y  = h(FX);

    %%% Calculate HF via affine regression from X to Y = h(f(X)) %%%
    
    yf = sum(repmat(weights,M,1).*Y,2);
    deltaX = X(:,2:end)-repmat(xa,1,2*N);
    deltaY = Y(:,2:end)-repmat(yf,1,2*N);
    HF = (deltaY/deltaX);

    %%% Calculate the a forecast state estimate, xf, and covariance, Pxx %%%
    
    xf = sum(repmat(weights,N,1).*FX,2);    
    deltaFX = FX(:,2:end)-repmat(xf,1,2*N);
    FPF = (deltaFX*weightMatrix*deltaFX');   
    Pxx = FPF+Q;   

    %%% Generate an updated unscented ensemble and recompute the observation, h %%%
    
    [FX,weights] = unscentedEnsemble(xf,Pxx,sqrt(N-1));
    deltaFX = FX(:,2:end)-repmat(xf,1,2*N); 
    Y = h(FX);
    
    %%% Calculate predicted observations and associated statistics %%%

    yf = sum(repmat(weights,M,1).*Y,2);
    deltaY = Y(:,2:end)-repmat(yf,1,2*N);
    Pb = deltaY*weightMatrix*deltaY';
    Py = Pb+R;
    Pxy = deltaFX*weightMatrix*deltaY';

    %%% Kalman Update %%%
    
    K = Pxy/Py;
    HT=(Pxx\Pxy);
    
    innov = (obs-yf);    
    xa = xf+K*innov;
    Pa = Pxx - K*Py*K'; 
 
    %%% Estimate the observation and system noise after the Kalman update %%%
    
    if (Q_param == 0)        
        Qest=0;Rest=0;
    else
        Gamma0 = innov*innov';
        Gamma0prev = prev_innov*prev_innov';
        Gamma1 = innov*prev_innov';

        Rest = Gamma0-Pb;
        if (Q_param == 1)       %%% Fit full Q matrix   
            M = (HF\Gamma1+K*Gamma0prev)/HT;
            Qest = M - prev_FPF;
        else                    %%% Fit a parameterized Q matrix
            C = (Gamma1+HF*K*Gamma0prev - HF*prev_FPF*HT);
            Qest = Qparamaterization(C(:),HF,HT,Q_param);
        end
    end


warning('on','MATLAB:rankDeficientMatrix');

 




    
    
    