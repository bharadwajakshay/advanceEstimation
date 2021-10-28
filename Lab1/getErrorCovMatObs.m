function [Cl] = getErrorCovMatObs(sigmaNot, input)
% Get error covariance matrix for the  measurements

% Get the no of satellites seen. The dimention of the covariance matrix
% will be a square matrix with no of satellites seen 
Cls = eye(size(input,1));

% calculate error variance for each satelllite
for i=1:size(input,1)
    Cl(i,i) = sigmaNot.^2/(sind(input(i,end)))^2;
end

end