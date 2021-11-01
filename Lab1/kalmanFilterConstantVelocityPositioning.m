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
stateVec = zeros(8);
stateVec(1:4) = initialEstimateKF;

P = eye(size(P,1));
Qvec = [1, 1, 1, 1]; 

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
       
       
end


