function [H, R, Po] = getHandRMat(currentdata,estimatedRecPos,sigmaNot)

    % Caluclate R and H for calculating Kalman Gain
    R = getErrorCovMatObs(sigmaNot, currentdata);
    
    Po = sqrt((currentdata(:,3)-estimatedRecPos(1)).^2 + ...
            (currentdata(:,4)-estimatedRecPos(2)).^2 +...
            (currentdata(:,5)-estimatedRecPos(3)).^2)+estimatedRecPos(4);
        
    H = [(estimatedRecPos(1) - currentdata(:,3))./Po, ...
            (estimatedRecPos(2) - currentdata(:,4))./Po, ...
            (estimatedRecPos(3) - currentdata(:,5))./Po, ones(size(Po))];
end

