function [ simLosses ] = loadSimLossFromExcel(absPathToExcel, ...
    MULTIPLICATION_FACTOR)
%LOADSIMLOSSFROMEXCEL A helper function to load the simulation path loss
%results from the input Excel file.
%
% Yaguang Zhang, Purdue, 10/31/2019

if ~exist('MULTIPLICATION_FACTOR', 'var')
    MULTIPLICATION_FACTOR = 1;
end

% Flag to control whether to load in the raw sim results or the calibrated
% one.
flagDataSource = 'raw';

simResultsTable = readtable(absPathToExcel);

switch lower(flagDataSource)
    case 'raw'
        
        refPeak = simResultsTable.simulatedPeak_dB_(1);
        simLosses = (refPeak - simResultsTable.simulatedPeak_dB_) ...
            .*MULTIPLICATION_FACTOR;
    case 'calibrated'
        simLosses = simResultsTable.simLoss;
    otherwise
end

end
% EOF