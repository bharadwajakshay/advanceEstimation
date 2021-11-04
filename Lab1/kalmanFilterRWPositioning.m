clc;
clear all;
close all;

%% Read the data from the file 
refPosFileName = "data/refpos.txt";
refPos = readmatrix(refPosFileName);

satPosFileName = "data/satpos_meas.txt";
satPos = readmatrix(satPosFileName);

%% Implementation of Kalman Filter
timeEpochs = unique(satPos(:,1));
sigmaNot = 1;
c = physconst('LightSpeed');

% Get least sq estimate for the 1st epoch 
% Set initial reciever position to zero. At the center of the earth
% Seed for the least squares estimate
estimatedRecPos = [0, 0, 0, 0]';
idx = find(satPos(:,1) == timeEpochs(1));
currentdata = satPos(idx,:);
Cl= getErrorCovMatObs(sigmaNot, currentdata);
initialEstimateKF = leastSquareEstimate(currentdata, Cl, estimatedRecPos, sigmaNot);

% Initialize the Kalman Filter 
% Contains [x,y,z,cdt]'
stateVec = initialEstimateKF;

% multiply dt with c to get the actual cdt
stateVec(4) = stateVec(4)*c;

% transitional Matrix
Phi = eye(4);
% Consider very high uncertainty in the position vector. 
P = eye(4) .* 10;

% Process noise 
Qvec = [1, 1, 1, 1]';


kalmanPos = zeros(size(refPos,1),3);
kalmanclockDriftErr = zeros(size(refPos,1),1);

kalmanPos(1,:) = stateVec(1:3);

h = waitbar(0,'Performing Kalman estimation');

for eachtimestamp = 2:size(timeEpochs,1)
    deltaT = timeEpochs(eachtimestamp,1) - timeEpochs(eachtimestamp-1,1);
    Q = diag(Qvec).*deltaT;
    
    % Start Kalman Filter
    % Pridiction
    stateVec = Phi * stateVec;
    P = Phi * P * Phi' + Q;
    
    % Calculate H and R matrix
    idx = find(satPos(:,1) == timeEpochs(eachtimestamp));
    currentdata = satPos(idx,:);
    
    % Get measurements
    Z = currentdata(:,6);
    
    [H, R, Po] = getHandRMat(currentdata, stateVec, sigmaNot);
    
    % Calculate Kalman Gain
    
    K = P*H'/(H*P*H'+R);
    
    % Update Step  
    stateVec = stateVec + K*(Z - Po);
    
    P = (eye(size(P,1))- (K*H))*P;  
    
    kalmanPos(eachtimestamp,:) = stateVec(1:3);
    
    kalmanclockDriftErr(eachtimestamp) = stateVec(4);

    waitbar(eachtimestamp/size(timeEpochs,1),h)
end

close(h)



%% Plot the Estimated positions
% Plot the reference positions
figure
plot3(refPos(:,2),refPos(:,3),refPos(:,4),'-');
grid on;
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('3D plot of the true trajectory vs Kalman filter estimate');
hold on;
plot3(kalmanPos(:,1),kalmanPos(:,2),kalmanPos(:,3),'-');
legend('True Reference Position', 'Kalman Filter Estimate')

figure;
truegeoposition = ecef2lla(refPos(:,2:4));
geoposition = ecef2lla(kalmanPos(:,1:3));
geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'b');
hold on;
geoplot(geoposition(:,1),geoposition(:,2),"Color",'r');
grid off;
legend('True Reference Position', 'Kalman Filter Estimate ')
title(['Plot of the latitude and longitude of true position agaist...' ...
    'the Kalman filter estimate'])

figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,kalmanclockDriftErr./c);
grid on;
title('Plot of the clock drift');
xlabel('Time epoch in sec');
ylabel('Time error in sec');

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



