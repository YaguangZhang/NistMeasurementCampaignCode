% GENERATEPDPPLOTS Generate the PDP overview .png plot for recorded signals
% and save some of them as .fig file for further inspectation.
%
% Yaguang Zhang, Purdue, 04/30/2018

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PdpInspection');

% Reuse results from loadMeasCampaignInfo.m.
ABS_PATH_TO_TX_INFO_LOGS_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

% Reuse results from evalPathLossesForContiTracks.m and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');

% Which plots to also save as .fig files, besides the default .png plots.
INDICES_MEAS_PLOTS_TO_SAVE_AS_FIG = [];

% Generate debug figures for noise elimination or not.
flagGenerateNoiseEliDebugFig = true;

%% Before Processing the Data

curFileName = mfilename;
fileNameHintRuler = [' ', repmat('-', 1, length(curFileName)+2), ' '];
disp(fileNameHintRuler)
disp(['  ', curFileName, '  '])
disp(fileNameHintRuler)

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% No need to show the figures since all of them are saved to the disk.
curDefaultFigureVisible = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');

%% Get Info for Measurement Data Files and Calibration Polynomials

disp(' ')
disp('    Loading results from: ')
disp('      - txInfoLogs.mat')
disp('      - contiPathLossesWithGpsInfo.mat')

assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');
assert(exist(ABS_PATH_TO_PATH_LOSSES, 'file')==2, ...
    'Couldn''t find contiPathLossesWithGpsInfo.mat! Please run PostProcessing/3_PathLossComputation/evalPathLossesForContiTracks.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')

% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Get 'contiPathLossesWithGpsInfo', 'contiOutFilesRelPathsUnderDataFolder'
% and 'contiOutFileIndicesReflection'.
load(ABS_PATH_TO_PATH_LOSSES);

% Parameters need for signal noise elimination.
Fs = F_S;

disp('    Done!')

%% Locate All Signal Recordings

disp(' ')
disp('    Locating signal recording files...')

absPathToMeasDataDir = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'Data', ...
    '20180331_NistFoliage');
allSigOutFilesMeas = locateFilteredOutFiles(absPathToMeasDataDir);

disp('    Done!')

%% Generate and Save PDP Plots

disp(' ')
disp('    Plotting PDPs...')

[timeNsInPlots, signalAmpInPlots, lowPassedSigsInPlot, ...
    plotFileNames, sampIndexRanges] ...
    = inspectPdps('Meas', allSigOutFilesMeas, ABS_PATH_TO_SAVE_PLOTS, ...
    INDICES_MEAS_PLOTS_TO_SAVE_AS_FIG, F_S, flagGenerateNoiseEliDebugFig);

disp('    Done!')

%% Estimate the Energy Ratio of LOS Signals

disp(' ')
disp('    Estimate the energy ratio of LOS signals for the PDPs...')

fullPathsToSaveEnRatioPlots = cellfun(@(p) fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [p, '_EnergyRatio.png']), plotFileNames, 'UniformOutput', false);
energyRatiosForLosSig ...
    = cellfun(@(ts, as, lps, p) ...
    estimateEnergyRatioInOnePdpForLosSig(ts, as, p, lps), ...
    timeNsInPlots, signalAmpInPlots, ...
    lowPassedSigsInPlot, fullPathsToSaveEnRatioPlots);

disp('    Done!')

%% Locate the PDPs on the Map

disp(' ')
disp('    Generating map plots for the PDPs...')

fullPathsToSaveMapPlots = cellfun(@(p) fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [p, '_LocOnMap.png']), plotFileNames, 'UniformOutput', false);
% Manually adjust the visible area.
axisToSet = [-105.2776 -105.2743 39.9892 39.9917];
numFigs = length(sampIndexRanges);
idxRelavantGpsSampVec = nan(numFigs,1);
for idxTrack = 1:numFigs
    [~, curIdxRelavantGpsSamp] = locateSigSegOnMap( ...
        contiPathLossesWithGpsInfo, contiPathLossesExtraInfo, ...
        idxTrack, sampIndexRanges{idxTrack}, ...
        fullPathsToSaveMapPlots{idxTrack}, axisToSet);
    assert(length(curIdxRelavantGpsSamp)==1, ...
        'More than one GPS samples were found for the PDP plot!');
    idxRelavantGpsSampVec(idxTrack) = curIdxRelavantGpsSamp;
end

save(fullfile(ABS_PATH_TO_SAVE_PLOTS, 'idxRelavantGpsSampVec'), ...
    'idxRelavantGpsSampVec');

disp('    Done!')

%% Set DefaultFigureVisible Back
set(0, 'DefaultFigureVisible', curDefaultFigureVisible);

% EOF