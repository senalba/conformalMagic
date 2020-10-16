function [cosI, RI] = cosVsR(qRadon, DRadon, idxPop, idxP, idxRP)
    
    unitDirection = cell(1,length(idxPop));
    meanRadonDistance = zeros(1,length(idxPop));
    
    for iP = idxP
        
        [~, ~, ~, unitDirection{iP}, meanRadonDistance(iP)] = ...
            centroidsMutants(qRadon, idxPop{iP}, idxRP);
    end
    
    
    [eigVectRP, ~, ~] = pca(qRadon(idxRP,:));
    
    
    for iP = 1:length(idxP)
        
        RI(iP) = meanRadonDistance(idxP(iP));
    end
    
    iPCA = 1;
    
    for iP = length(idxP):-1:1
        
        cosI(iP) = abs(dot( unitDirection{idxP(iP)}, eigVectRP(:,iPCA) ));
    end    
    