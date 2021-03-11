function forecast = NonParametricDynamics(state,trainingData,kNN,linearModel,useWeights,forecastHorizon)
% state             N-by-E matrix of E state vectors in R^N to forecast
% trainingData      N-by-T matrix of T state vectors in R^N used for forecasting
% kNN               Number of nearest neighbors for the forecast
% linearModel       1 - Use a locally linear model, 0 - use a local average
% useWeights        1 - weight neighbors based on distance, 0 - no weights
% forecastHorizon   Number of steps to forecast, use 1 for filtering
    
    if (nargin < 6) forecastHorizon = 1;    end
    if (nargin<5)   useWeights = 1;         end 
    if (nargin<4)   linearModel = 0;        end  
    if (nargin<3)   kNN = 40;               end 
    
    %%% Find the nearest neighbors
    [ds,inds] = pdist2(trainingData(:,1:end-forecastHorizon)',state','euclidean','smallest',kNN);

    if useWeights  
        sigma = (1/2)*mean(ds);
        weights = exp(-(ds./repmat(sigma,size(ds,1),1)));
        weights = weights./repmat(sum(weights),size(weights,1),1);
    else
        weights = ones(kNN,size(inds,2))/kNN;
    end
 
    forecast = zeros(size(state));
    
    if linearModel
        for i = 1:size(inds,2)
            Fnn = trainingData(:,inds(:,i)+forecastHorizon);%%% targets
            muFnn = Fnn*weights(:,i);                       %%% target mean
            nn = trainingData(:,inds(:,i));                 %%% source points
            munn = nn*weights(:,i);                         %%% source mean
            llModel = (Fnn-repmat(muFnn,1,size(Fnn,2)))/(nn - repmat(munn,1,size(nn,2)));
            forecast(:,i) = llModel*(state(:,i) - munn) + muFnn; %%% linear model
        end     
    else
        for i = 1:size(inds,2)    %%% weighted average of neighbor forecasts
            forecast(:,i) = trainingData(:,inds(:,i)+forecastHorizon)*weights(:,i);
        end
    end
    
end


