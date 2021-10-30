clc;
clear all;
close all;

%% Read the data from the file 
refPosFileName = "data/refpos.txt";
refPos = readmatrix(refPosFileName);

satPosFileName = "data/satpos_meas.txt";
satPos = readmatrix(satPosFileName);

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
% Plot the reference positions
plot3(refPos(:,2),refPos(:,3),refPos(:,4),'-');
grid on;
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('True position');
hold on;
plot3(leaseSqPos(:,1),leaseSqPos(:,2),leaseSqPos(:,3),'-');
legend('True Reference Position', 'Least Squares Estimate')

figure;
truegeoposition = ecef2lla(refPos(:,2:4));
geoposition = ecef2lla(leaseSqPos(:,1:3));
geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'k');
hold on;
geoplot(geoposition(:,1),geoposition(:,2),"Color",'g');
grid off;
legend('True Reference Position', 'Least Squares Estimate')
title(['Plot of the latitude and longitude of true position agaist...' ...
    'the least square estimate'])

figure;
timeEpochsPlot = timeEpochs-timeEpochs(1);
plot(timeEpochsPlot,clockDrift);
grid on;

%% Implementation of Kalman Filter

% Get least sq estimate for the 1st epoch 

% Set initial reciever position to zero. At the center of the earth
% Seed for the least squares estimate
estimatedRecPos = [0, 0, 0, 0]';
idx = find(satPos(:,1)==timeEpochs(1));
currentdata = satPos(idx,:);
Cl= getErrorCovMatObs(sigmaNot, currentdata);
initialEstimateKF = leastSquareEstimate(currentdata, Cl, estimatedRecPos, sigmaNot);

% Initialize the Kalman Filter 
% Contains [x,y,z,cdt]'
stateVec = initialEstimateKF;

% transitional Matrix
Phi = eye(4);
% Consider very high uncertainty in the position vector. 
P = eye(4) .* 10;

% Process noise 
Qvec = [1, 1, 1, 1]';
c = physconst('LightSpeed');

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
    
    [H, R] = getHandRMat();
    
        
    % Calculate Kalman Gain
    
    K = P*H'*inv(H*P*H'+R);
    
    % Update Step
    
    % Get measurements
    Z = currentdata(:,6);
    
    stateVec = stateVec + K*(Z - (H*stateVec));
    
    P = (eye(size(P,1))- (K*H))*P;  
    
    kalmanPos(eachtimestamp,:) = stateVec(1:3);

    waitbar(eachtimestamp/size(timeEpochs,1),h)
end



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

% figure;
% truegeoposition = ecef2lla(refPos(:,2:4));
% geoposition = ecef2lla(leaseSqPos(:,1:3));
% geoplot(truegeoposition(:,1),truegeoposition(:,2),"Color",'k');
% hold on;
% geoplot(geoposition(:,1),geoposition(:,2),"Color",'g');
% grid off;
% legend('True Reference Position', 'Least Squares Estimate')
% title(['Plot of the latitude and longitude of true position agaist...' ...
%     'the least square estimate'])
% 
% figure;
% timeEpochsPlot = timeEpochs-timeEpochs(1);
% plot(timeEpochsPlot,clockDrift);
% grid on;






  




