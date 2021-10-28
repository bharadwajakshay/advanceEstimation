clc;
clear all;
close all;

%% Read the data from the file 
refPosFileName = "data/refpos.txt";
refPos = readmatrix(refPosFileName);

satPosFileName = "data/satpos_meas.txt";
satPos = readmatrix(satPosFileName);

%% Plot the reference positions
plot3(refPos(:,2),refPos(:,3),refPos(:,4),'-');
grid on;
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('True position');

%% implementation of LS Method
timeEpochs = refPos(:,1);

sigmaNot = 1;

c = physconst('LightSpeed');

deltaPs = ones(1,3)*100;

% Set initial reciever position to zero. At the center of the earth
estimatedRecPos = [0, 0, 0, 0]';

leaseSqPos = zeros(size(refPos,1),3);
clockDrift = zeros(size(refPos,1),1);

h = waitbar(0,'Performing least squares estimate');

for eachtimestamp = 1:size(timeEpochs,1)

    idx = find(satPos(:,1)==timeEpochs(eachtimestamp));
    
    % Data needed to be used for estimation for current epoch
    currentdata = satPos(idx,:);

    % Calculate the Covariance matrix for measurement
    Cl= getErrorCovMatObs(sigmaNot, currentdata);

    stop = false;

    count = 0;
    prevMean = 100;
    while ~stop
        % Step 1: caluclate the pseudo ranges for each satellite
        % measured pseudoranges
        Ps = currentdata(:,6);

        % Step 2: Calculate the pseudo ranges based on the estimated value
        Po = sqrt((currentdata(:,3)-estimatedRecPos(1)).^2 + ...
            (currentdata(:,4)-estimatedRecPos(2)).^2 +...
            (currentdata(:,5)-estimatedRecPos(3)).^2);

        % Step 2.a: Calculate the delta Pseudo ranges
        deltaP = Ps - Po;

        % Step 3: Calculate the design matrix
        A = [(estimatedRecPos(1) - currentdata(:,3))./Po, ...
            (estimatedRecPos(2) - currentdata(:,4))./Po, ...
            (estimatedRecPos(3) - currentdata(:,5))./Po, repmat(c,(size(Po)))];

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

         if abs(prevMean - mean(deltaPs)) < 0.01
             stop = true;
         else
             prevMean = mean(deltaPs);
         end

        count = count + 1;
    end

    
    leaseSqPos(eachtimestamp,:) = estimatedRecPos(1:3)';
    clockDrift(eachtimestamp) = estimatedRecPos(4);

    waitbar(eachtimestamp/size(timeEpochs,1),h)

end
close(h);

%% Plot the Estimated positions
hold on;
plot3(leaseSqPos(:,1),leaseSqPos(:,2),leaseSqPos(:,3),'-');
legend('True Reference Position', 'Least Squares Estimate')

figure;
timeEpochs = timeEpochs-timeEpochs(1);
plot(timeEpochs,clockDrift);
grid on;

%% Implementation of Kalman Filter

% Get least sq estimate for the 1st epoch 

disp('Breakpoint')




