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
    
    % Perform least squares estimate
    estimatedRecPos = leastSquareEstimate(currentdata, Cl, estimatedRecPos,sigmaNot);

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




