function [truth,obs,p] = GenerateL96(T,N,dt,Q,R)
    
    truth = zeros(N,T);     %%% noise free 'truth'
    obs = zeros(N,T);       %%% noisy observations
    
    transient = 1000;
    substeps = ceil(dt/.05);
    h = dt/substeps;
    state = rand(N,1);
    
    %%% Parameters in the L96 vector field
    a = 1*ones(N,1);
    F = 8*ones(N,1);
    p = [a F];
    
    %%% Compute the matrix square roots of Q and R, used to generate
    %%% appropriately correlated noise samples below
    [U,S,V] = svd(Q);
    rootQ = U*diag(sqrt(diag(S)))*U';
    [U,S,V] = svd(R);
    rootR = U*diag(sqrt(diag(S)))*U';
    
    %%% Run an initial transient to get onto the attractor
    for i = 1:transient
        state = L96Dynamics(state,dt,p);
    end
    
    for i = 1:T
        for k = 1:substeps
            %%% RK4
            k1=h*LorenzVectorField(state,N,p);
            k2=h*LorenzVectorField(state+k1/2,N,p);
            k3=h*LorenzVectorField(state+k2/2,N,p);
            k4=h*LorenzVectorField(state+k3,N,p);
            state=state+k1/6+k2/3+k3/3+k4/6;
            
            %%% Stochastic forcing with coviariance Q
            state = state+sqrt(h)*rootQ*randn(N,1);
        end
        %%% Clean 'true' state
        truth(:,i) = state;
        %%% Noisy observed state with noise covariance R
        obs(:,i) = state+rootR*randn(N,1);
    end   
end



