function [ powerShiftsForCali ] ...
    = genCalibrationFct( lsLinesPolysInv, rxGains, gains )
%GENCALIBRATIONFCT Compute the power shifts required for converting
%calculated power(s) to compute measured power(s), according to the
%calibration line(s) for the input gain(s).
%
%  Inputs:
%    - lsLinesPolysInv
%      The inverse version of lsLinesPolys saved by calibrateRx.m. Note
%      that the lsLinesPolys are essentially for calculated power vs.
%      measured power. Here, we need measured power vs. calculated power.
%    - rxGains
%      Thr RX gains correspond to the polynormials in lsLinesPolysInv.
%    - gains
%      A scalar / a vector. Specify which calibration line is (or lines
%      are) needed.
%  Output:
%    - powerShiftsForCali
%      Essentially b in the linear polynomial:
%          measuredPowerInDb = caluculatedPowerInDb + b,
%      since we will force it to use 1 as its slope.
%
% Update 04/10/2018: Minor fix for NIST dataset (with only one calibration
% line available).
%
% Yaguang Zhang, Purdue, 09/14/2017

switch length(rxGains)
    case 1
        if ~all(gains==rxGains)
            warning(['    ', num2str(sum(gains~=rxGains)), ...
                ' of samples were collected with gain not clibratable!'])
        end
        
        powerShiftsForCali = ones(length(gains),1).*lsLinesPolysInv{1}(2);
    case 2
        % lsLinesPolysInv{1} is for the first calibration data set.
        powerShiftsForCali = ((gains-rxGains(2)).*lsLinesPolysInv{1}(2) ...
            +(rxGains(1)-gains).*lsLinesPolysInv{2}(2)) ...
            ./(rxGains(1)-rxGains(2));
    otherwise
        error('Unexpected number of calibration lines found in lsLinesPolysInv!')
end

end
% EOF