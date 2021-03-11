function obs = L96Obs(state,M)
%%% M is the dimension of the observation, N=M is full obs
    N = size(state,1);
    spacing = ceil(N/M);
    
    obs = state(1:spacing:N,:);

end

