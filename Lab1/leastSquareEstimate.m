function [estimatedRecPos,deltaPs] = leastSquareEstimate(currentdata, Cl, estimatedRecPos, sigmaNot)

    c = physconst('LightSpeed');
    stop = false;
    count = 0;
    prevMean = 100;
    
    deltaPs = ones(5).*1000000;

    while ~stop
        % Step 1: caluclate the pseudo ranges for each satellite
        % measured pseudoranges
        Ps = currentdata(:,6);

        % Step 2: Calculate the pseudo ranges based on the estimated value
        Po = sqrt((currentdata(:,3)-estimatedRecPos(1)).^2 + ...
            (currentdata(:,4)-estimatedRecPos(2)).^2 +...
            (currentdata(:,5)-estimatedRecPos(3)).^2);

        % Step 2.a: Calculate the delta Pseudo ranges
        deltaP = Ps - (Po+(estimatedRecPos(4)*c));

        % Step 3: Calculate the design matrix
        A = [(estimatedRecPos(1) - currentdata(:,3))./Po, ...
            (estimatedRecPos(2) - currentdata(:,4))./Po, ...
            (estimatedRecPos(3) - currentdata(:,5))./Po, repmat(c,size(Po))];

        % Step 3.a: Calculate the weight matrix
        if (det(Cl)~=0)
            P = sigmaNot^2 * inv(Cl);
        else
            disp('The Cl matrix is singular. So calulating the pseudoInvers')
            P = sigmaNot^2 * pinv(Cl);
        end

        % Step 4: Perform lease sq estimation
        deltaX = inv(A'*P*A)*A'*P*deltaP;
        % Step 5: Check if the termination condition is met

        estimatedRecPos = estimatedRecPos + deltaX;

        % Euclidean distance between Estimated position and refrence
        % position

        % Step 2: ReCalculate DeltaP
        Po = sqrt((currentdata(:,3)-estimatedRecPos(1)).^2 + ...
            (currentdata(:,4)-estimatedRecPos(2)).^2 +...
            (currentdata(:,5)-estimatedRecPos(3)).^2);

        % Step 2.a: Calculate the delta Pseudo ranges
        deltaP = Ps - Po;

        errorCriteria = deltaP'*Cl*deltaP;

%         euclideanDist = sqrt((refPos(eachtimestamp,2)- estimatedRecPos(1))^2 + ...
%             (refPos(eachtimestamp,3)- estimatedRecPos(2))^2 + ...
%             (refPos(eachtimestamp,4)- estimatedRecPos(3))^2);

%        deltaPs(mod(count,3)+1) = euclideanDist;
        
         deltaPs(mod(count,3)+1) = errorCriteria;

         if abs(prevMean - mean(deltaPs)) < 0.00001
             stop = true;
         else
             prevMean = mean(deltaPs);
         end

        count = count + 1;
        
    end
end