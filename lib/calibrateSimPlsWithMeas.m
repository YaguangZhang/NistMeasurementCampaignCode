function [calibratedSimPls, shift, multiFactor] ...
    = calibrateSimPlsWithMeas(simPls, measPls, method)
%CALIBRATESIMPLSWITHMEAS Calibrate the simulation path loss results
%according to the corresponding measurements.
%
% Inputs:
%   - simPls, measPls
%     The simulation path losses and the reference measurement results.
%
% Outputs:
%   - calibratedSimPls
%     The calibrated simulation results.
%   - shift, multiFactor
%     The resulting parameters which are used via:
%           calibratedSimPls = simPls.*multiFactor+shift .
%   - method
%     Optional flag ('both' or 'shiftOnly') to control whether it is only
%     to use multiFactor~=1.
%
% Yaguang Zhang, Purdue, 10/31/2019

if ~exist('method', 'var')
    method = 'both';
end

expectedNumOfSamps = length(measPls);
assert(length(simPls)==expectedNumOfSamps, ...
    'Unexpected number of simulation results!');

mseFct = @(shift) sum((simPls+shift-measPls).^2)/expectedNumOfSamps;

switch lower(method)
    case 'both'
        fctToFit = @(paraToFit, input) input.*abs(paraToFit(1)) ...
            + paraToFit(2);
        [~, shiftStart, ~] = calibrateSimPlsWithMeas( ...
            simPls, measPls, 'ShiftOnly');
        
        startingPt = [30 shiftStart];
        
        fittedRes = nlinfit(simPls, measPls, fctToFit, startingPt);
        
        multiFactor = abs(fittedRes(1));
        shift = fittedRes(2);
    case 'shiftonly'
        minShift = min(measPls)-max(simPls);
        shift = fminsearch(mseFct, minShift);
        multiFactor = 1;
    otherwise
        error(['Unknown calibration method ', method, '!']);
end

calibratedSimPls = simPls.*multiFactor + shift;

end
% EOF