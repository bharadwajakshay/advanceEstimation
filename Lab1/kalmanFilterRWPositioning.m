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
    stateVec = stateVec + K*(Z - (H*stateVec));
    
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
title('True position');
hold on;
plot3(kalmanPos(:,1),kalmanPos(:,2),kalmanPos(:,3),'-');
legend('True Reference Position', 'Least Squares Estimate')

figure;
truegeoposition = ecef2lla(refPos(:,2:4));
geoposition = ecef2lla(kalmanPos(:,1:3));
geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'k');
hold on;
geoplot(geoposition(:,1),geoposition(:,2),"Color",'g');
grid off;
legend('True Reference Position', 'Least Squares Estimate')
title(['Plot of the latitude and longitude of true position agaist...' ...
    'the least square estimate'])

figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,kalmanclockDriftErr./c);
grid on;






  




