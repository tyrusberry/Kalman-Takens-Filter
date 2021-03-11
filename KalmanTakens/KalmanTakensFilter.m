function [filteredTimeSeries,var,Q,R] = KalmanTakensFilter(timeSeries,delays,kNN,Q_param,stationarity)
% timeSeries    a N-by-T time series to be filtered
% delays        number of delays to use in the Takens embeddings
% kNN           number of nearest neighbors to use in forecast
% Q_param       0 - none, 1 - full, 2 - multiple of identity, 3 - diagonal           
% stationarity  widow length for adaptively averaging Q and R

    if (nargin<5) stationarity = 2000; end
    if (nargin<4) Q_param = 1; end

    N=size(timeSeries,1);
    trainingData = TakensEmbedding(timeSeries,delays);

    f=@(x) NonParametricDynamics(x,trainingData,kNN);   % no model dynamics
    h=@(x) NonParametricObs(x,delays);                  % no model obs function
    state = trainingData(:,1);                          % initial state

    % run the filter
    [filteredTimeSeries,var,Q,R] = AUKF(state,timeSeries,f,h,Q_param,stationarity);    
    filteredTimeSeries = filteredTimeSeries(1:N,:);

end

