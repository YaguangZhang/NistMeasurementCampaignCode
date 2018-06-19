% POSTPROCESSMEASDATA Post process the NIST measurement dataset.
%
% This is a holder script listing all the steps required to properly
% process the measurement dataset. Please comment & uncomment commands as
% it is needed, depending on which results require updates.
%
% Yaguang Zhang, Purdue, 06/08/2018

clear; clc; close all;

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(pwd)));
setPath;

%% 1_Calibration: Calibrate the Gnu Radio RX
addpath(fullfile(pwd, '1_Calibration'));
calibrateRx;

%% 2_PlotGpsSampsOnMap: Plot Overview Maps
addpath(fullfile(pwd, '2_PlotGpsSampsOnMap'));
plotAllGpsSampsOnMap;

%% 3_PathLossComputation: Compute & Plot the Path Losses

% This will generate: 'TX_POWER_DBM', 'TX_HEIGHT_FEET', 'TX_HEIGHT_M',
% 'F_S', 'TX_LAT', 'TX_LON', 'TX_INFO_LOGS', 'TX_INFO_LOGS_ABS_PAR_DIRS', etc.
addpath(fullfile(pwd, '3_PathLossComputation'));
loadMeasCampaignInfo;

% Set this to be true to get rid of data too close to the TX.
flagTailorData = true;
% This will actually carry out the path loss computation.
if flagTailorData
    addpath(fullfile(pwd, '3_1_PathLossComputationTailored'));
    evalPathLossesForContiTracksTailored;
else
    evalPathLossesForContiTracks;
end

%% 4_FoliageAttenuationEstimation: Estimate Excess Path Losses
addpath(fullfile(pwd, '4_FoliageAttenuationEstimation'));
estimateFoliageAttenuation;

%% 5_PdpInspection: Plot One PDP per Track
addpath(fullfile(pwd, '5_PdpInspection'));
generatePdpPlots;

%% 6_DebugContiPathLosses: Plot More PDPs
debugContiPathLosses;

%% 7_ManuallyLocateTrees: Manually Label Trees on Map

% Set this to be true if it is necessary to re-label the trees.
flagRelabelTrees = false;
if flagRelabelTrees
    addpath(fullfile(pwd, '7_ManuallyLocateTrees'));
    manuallyLocateTrees;
end

%% 8_FoliageAttenuationEstimation_ManualTreeLocs
addpath(fullfile(pwd, '8_FoliageAttenuationEstimation_ManualTreeLocs'));
estimateFoliageAttenuationWithManualTreeLocs;

% EOF