clc;
clear all;
close all;

%% Read the data from the file 
refPosFileName = "data/refpos.txt";
refPos = readmatrix(refPosFileName);

satPosFileName = "data/satpos_meas.txt";
satPos = readmatrix(satPosFileName);

%% implementation of LS Method
timeEpochs = unique(satPos(:,1));

sigmaNot = 1;

deltaPs = ones(1,3)*100;

% Set initial reciever position to zero. At the center of the earth
estimatedRecPos = [0, 0, 0, 0]';

leaseSqPos = zeros(size(refPos,1),3);
clockDriftErr = zeros(size(refPos,1),1);

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
    clockDriftErr(eachtimestamp) = estimatedRecPos(4);

    waitbar(eachtimestamp/size(timeEpochs,1),h)

end
close(h);

%% Plot the Estimated positions
% Plot the reference positions
plot3(refPos(:,2),refPos(:,3),refPos(:,4),'-');
grid on;
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('3D plot of the true trajectory vs least sq estimate');
hold on;
plot3(leaseSqPos(:,1),leaseSqPos(:,2),leaseSqPos(:,3),'-');
legend('True Reference Position', 'Least Squares Estimate')

figure;
truegeoposition = ecef2lla(refPos(:,2:4));
geoposition = ecef2lla(leaseSqPos(:,1:3));
geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'b');
hold on;
geoplot(geoposition(:,1),geoposition(:,2),"Color",'r');
grid off;
legend('True Reference Position', 'Least Squares Estimate')
title(['Plot of the latitude and longitude of true position agaist...' ...
    'the least square estimate'])

figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,clockDriftErr);
title('Plot of the clock drift');
xlabel('Time epoch in sec');
ylabel('Time error in sec');
grid on;

%%  Error Analysis

truePositionENU = lla2enu(truegeoposition,truegeoposition(1,:),'ellipsoid');
estPos = lla2enu(geoposition,truegeoposition(1,:),'ellipsoid');
euclideanDistance = sqrt((truePositionENU(:,1)-estPos(:,1)).^2+...
                         (truePositionENU(:,2)-estPos(:,2)).^2+...
                         (truePositionENU(:,3)-estPos(:,3)).^2);
meanPosErr = mean(euclideanDistance);
stdPosErr = std(euclideanDistance);

disp('Mean Position error is: ')
disp(meanPosErr)
disp('Standar deviation in Position error is: ')
disp(stdPosErr)

figure;
plot(timeEpochsPlot,euclideanDistance);
title('Position error at each time epoch')
xlabel('Time in secs')
ylabel('Error in m')
grid on;

figure;
subplot(3,1,1);
plot(timeEpochsPlot,truePositionENU(:,1)-estPos(:,1))
subtitle('Error in East direction');
ylabel('Error in m')
grid on;
title('Position error in individual axes in ENU frame')

subplot(3,1,2);
plot(timeEpochsPlot,truePositionENU(:,2)-estPos(:,2))
subtitle('Error in North direction');
ylabel('Error in m')
grid on;

subplot(3,1,3);
plot(timeEpochsPlot,truePositionENU(:,3)-estPos(:,3))
subtitle('Error in Up direction');
ylabel('Error in m')
grid on;

