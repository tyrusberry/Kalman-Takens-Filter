function b=rema(a,N)
    %%% Lorenz-96 takes place on a circle of N nodes, this computes the
    %%% appropriate node for inputs outside of 1,...,N
    b=mod(a-1,N)+1;
end

