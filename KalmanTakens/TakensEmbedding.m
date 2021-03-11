function delaySeries = TakensEmbedding(timeSeries,delays)
    
    delaySeries = [];
    
    for i = delays:-1:1
        
        delaySeries = [delaySeries; timeSeries(:,i:end-(delays-i))];
        
    end
   
    





