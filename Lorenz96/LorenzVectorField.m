function ydot = LorenzVectorField(state,N,p)
    %%% Lorenz-96 vector field
    a = p(:,1);
    F = p(:,2);
    ydot=zeros(size(state));
    for i=1:N
        ydot(i,:)=a(i)*state(rema(i+1,N),:).*state(rema(i-1,N),:)-state(rema(i-2,N),:).*state(rema(i-1,N),:)-state(rema(i,N),:)+F(i);
    end
end


