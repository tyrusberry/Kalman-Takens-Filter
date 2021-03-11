function obs = NonParametricObs(state,delays)

    obs = state(1:delays:end,:); 
    