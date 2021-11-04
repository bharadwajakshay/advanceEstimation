clc;
clear all;
close all;

%% Read the data from the file 
refPosFileName = "data/refpos.txt";
refPos = readmatrix(refPosFileName);

satPosFileName = "data/satpos_meas.txt";
satPos = readmatrix(satPosFileName);


%% Implementation of Kalman Filter with constant velocity model

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
% multiply dt with c to get the actual cdt
initialEstimateKF(4) = initialEstimateKF(4)*c;

% Initialize the Kalman Filter 
% Contains [x,y,z,cdt,x.,y.,z.,cdt.]'
stateVec = zeros(8,1);
stateVec(1:4) = initialEstimateKF;

P = eye(size(stateVec,1));
Qvec = [1, 1, 1, 1]; 

kalmanPos = zeros(size(refPos,1),3);
kalmanclockDriftErr = zeros(size(refPos,1),1);
kalmanPosVelocity = zeros(size(refPos,1),3);
kalmanclockDriftErrVelocity = zeros(size(refPos,1),1);

kalmanPos(1,:) = stateVec(1:3);
kalmanPosVelocity(1,:) = stateVec(5:7);


h = waitbar(0,'Performing Kalman estimation');

for eachtimestamp = 2:size(timeEpochs,1)
    deltaT = timeEpochs(eachtimestamp,1) - timeEpochs(eachtimestamp-1,1);
    
    phi = [1 0 0 0 deltaT 0 0 0;
           0 1 0 0 0 deltaT 0 0;
           0 0 1 0 0 0 deltaT 0;
           0 0 0 1 0 0 0 deltaT;
           0 0 0 0 1 0 0 0;
           0 0 0 0 0 1 0 0;
           0 0 0 0 0 0 1 0;
           0 0 0 0 0 0 0 1];
    Q = [(Qvec(1)*deltaT^3/3) 0 0 0 (Qvec(1)*deltaT^2/2) 0 0 0;
         0 (Qvec(2)*deltaT^3/3) 0 0 0 (Qvec(2)*deltaT^2/2) 0 0;
         0 0 (Qvec(3)*deltaT^3/3) 0 0 0 (Qvec(3)*deltaT^2/2) 0;
         0 0 0 (Qvec(4)*deltaT^3/3) 0 0 0 (Qvec(4)*deltaT^2/2);
         (Qvec(1)*deltaT^2/2) 0 0 0 Qvec(1)*deltaT 0 0 0;
         0 (Qvec(2)*deltaT^2/2) 0 0 0 Qvec(2)*deltaT 0 0;
         0 0 (Qvec(3)*deltaT^2/2) 0 0 0 Qvec(3)*deltaT 0;
         0 0 0 (Qvec(4)*deltaT^2/2) 0 0 0 Qvec(4)*deltaT];
    
    stateVec = phi * stateVec;
    P = phi * P * phi' + Q;
    
    % Calculate H and R matrix
    idx = find(satPos(:,1) == timeEpochs(eachtimestamp));
    currentdata = satPos(idx,:);
    
    % Get measurements
    Z = currentdata(:,6);

    [H, R, Po] = getHandRMatCV(currentdata, stateVec, sigmaNot,deltaT);

    % Calculate Kalman Gain
    K = P*H'/(H*P*H'+R);

    % Update Step  
    stateVec = stateVec + K*(Z - Po);
    
    P = (eye(size(P,1))- (K*H))*P;

    kalmanPos(eachtimestamp,:) = stateVec(1:3);
    kalmanPosVelocity(eachtimestamp,:) = stateVec(5:7);
    kalmanclockDriftErr(eachtimestamp) = stateVec(4);
    kalmanclockDriftErrVelocity(eachtimestamp) = stateVec(8);




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
title('3D plot of the true trajectory vs Kalman filter estimate with constant velocity model');
hold on;
plot3(kalmanPos(:,1),kalmanPos(:,2),kalmanPos(:,3),'-');
legend('True Reference Position', 'Kalman Filter Estimate (Const Vel)')

figure;
truegeoposition = ecef2lla(refPos(:,2:4));
geoposition = ecef2lla(kalmanPos(:,1:3));
geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'b');
hold on;
geoplot(geoposition(:,1),geoposition(:,2),"Color",'r');
grid off;
legend('True Reference Position', 'Least Squares Estimate')
title(['Plot of the latitude and longitude of true position agaist...' ...
    ' Kalman filter estimate with constant velocity model'])

figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,kalmanclockDriftErr./c);
grid on;
title('Plot of the clock drift');
xlabel('Time epoch in sec');
ylabel('Time error in sec');

figure;
subplot(3,1,1);
plot(timeEpochsPlot,kalmanPosVelocity(:,1));
grid on;
subtitle('Velocity along X')
xlabel('Time in sec')
ylabel('Velocity in m/s')


subplot(3,1,2);
plot(timeEpochsPlot,kalmanPosVelocity(:,2));
grid on;
subtitle('Velocity along Y')
xlabel('Time in sec')
ylabel('Velocity in m/s')


subplot(3,1,3);
plot(timeEpochsPlot,kalmanPosVelocity(:,3));
grid on;
subtitle('Velocity along Z')
xlabel('Time in sec')
ylabel('Velocity in m/s')

title('Velocity plots')


figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,kalmanclockDriftErrVelocity./c);
grid on; 
title('Plot of the clock drift');
xlabel('Time epoch in sec');
ylabel('Time drift in sec');


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








