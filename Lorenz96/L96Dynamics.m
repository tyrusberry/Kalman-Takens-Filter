function state = L96Dynamics(state,dt,p)

    N = size(state,1);
    substeps = ceil(dt/.05);
    h = dt/substeps;
    
    for k = 1:substeps
        %%% RK4
        k1=h*LorenzVectorField(state,N,p);
        k2=h*LorenzVectorField(state+k1/2,N,p);
        k3=h*LorenzVectorField(state+k2/2,N,p);
        k4=h*LorenzVectorField(state+k3,N,p);
        state=state+k1/6+k2/3+k3/3+k4/6;
    end


end

