function [ pathLossInDb, absPathOutFile ] ...
    = computePathLossForOutFileDir(curOutFileDir, rxGain, ...
    noiseEliminationFct, powerShiftsForCali)
%COMPUTEPATHLOSSFOROUTFILEDIR Load the Gnu Radio samples stored in the .out
%file specified by the input dir struct outFileDir, and compute the path
%loss for it.
%
% We will consider both the TX calibration and the antenna normalization.
%
% Inputs:
%   - curOutFileDir
%     A dir struct specifying which .out file will be processed. Note:
%     curOutFileDir should at least has the field name, which contains the
%     full absolute path for the .out file to be processed, or only the
%     file name; For the second case, another field folder is required for
%     the full abosolute path of the parent folder.
%   - Other inputs
%     Please refer to computePathLossForCurSignal.m for more detail.
%
% Yaguang Zhang, Purdue, 09/26/2017

%% Parameters

% TX power (after the upconverter) in dBm.
try
    txPower = evalin('base', 'txPower');
catch
    warning('TX power (after the upconverter) txPower not found in the base workspace.')
    warning('Will use the default value -23 dBm.')
    txPower  = -23;
end

%% Load Data

% Generate the full path for the .out file.
if isfield(curOutFileDir, 'folder')
    [~, outFileName] = fileparts(curOutFileDir.name);
    absPathOutFile = fullfile(curOutFileDir.folder, [outFileName, '.out']);
else
    absPathOutFile = curOutFileDir.name;
end
curSignal = read_complex_binary(absPathOutFile);

%% Compute the Path Loss

pathLossInDb = computePathLossForCurSignal(curSignal, txPower, ...
    rxGain, noiseEliminationFct, powerShiftsForCali);

end
% EOF