% DEBUGCONTIPATHLOSSES Debug the continous path loss computation process by
% showing PDP's and antenna gain evaluation plots within a distance range.
%
% Yaguang Zhang, Purdue, 05/02/2018

clear; clc; close all;

%% Configurations

warning('on'); dbstop if error;

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% We will take advantage of the scripts for inspecting PDPs.
addpath(fullfile(pwd, '5_PdpInspection'));

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'ContiPathLossesDebugging');

% Reuse results from evalPathLossesForContiTracks.m and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');
ABS_PATH_TO_UTM_INFO= fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'FoliageAttenuationEstimation', ...
    'utmInfoForPathLossesAndTrees.mat');

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

%% Get Info for Measurement Data Files and Calibration Polynomials

disp(' ')
disp('    Loading results from: ')
disp('      - contiPathLossesWithGpsInfo.mat')
disp('      - txInfoLogs.mat')
disp('      - utmInfoForPathLossesAndTrees.mat')

assert(exist(ABS_PATH_TO_PATH_LOSSES, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/3_PathLossComputation/evalPathLossesForContiTracks.m first.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');
assert(exist(ABS_PATH_TO_UTM_INFO, 'file')==2, ...
    'Couldn''t find utmInfoForPathLossesAndTrees.mat! Please run PostProcessing/4_FoliageAttenuationEstimation/estimateFoliageAttenuation.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'contiPathLossesWithGpsInfo', 'contiOutFilesRelPathsUnderDataFolder'
% and 'contiOutFileIndicesReflection'.
load(ABS_PATH_TO_PATH_LOSSES);
% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Get location records in the UTM system, e.g. pathLossUtmXYHs,
% pathLossUtmZones, xTx, yTx, txUtmZone, treeLocations, treeUtmXYHs, and
% treeUtmZones.
load(ABS_PATH_TO_UTM_INFO);

% Parameters need for signal noise elimination.
Fs = F_S;
NUM_SIGMA_FOR_THRESHOLD = 3;

disp('    Done!')

%% Inspect PDP's within Specified Diantance Ranges

disp('    Generating PDP plots for signal recording segments ...');

% The first range contains everything. Other ranges are duplicated for more
% convenient analysis.
distRangesToInsp = {[50, 500], [0, 50], [175, 180]}; % In meter.
numDistRangesToInsp = length(distRangesToInsp);

losPeakEnergyRs = cell(numDistRangesToInsp, 1);
for idxRange = 1:numDistRangesToInsp
    curRange = distRangesToInsp{idxRange};
    curNumTracks = length(contiPathLossesWithGpsInfo);
    
    losPeakEnergyRsTracks = cell(curNumTracks, 1);
    disp(['        Range: ', ...
        num2str(idxRange), '/', num2str(length(distRangesToInsp))]);
    
    for idxTrack = 1:curNumTracks
        [numCurSegs, ~] = size(contiPathLossesWithGpsInfo{idxTrack});
        
        losPeakEnergyRsSegs = nan(numCurSegs, 1);
        disp(['            Track: ', ...
            num2str(idxTrack), '/', ...
            num2str(length(contiPathLossesWithGpsInfo))]);        
        
        % Take advantage of parallel computing for speed.
        curPathLossUtmXYHs = pathLossUtmXYHs{idxTrack};
        curContiOutFilesRelPathsUnderDataFolder ...
            = contiOutFilesRelPathsUnderDataFolder{idxTrack}{1};
        parfor idxSeg = 1:numCurSegs
            % Make Fs and USRP_NOISE_FLOOR_V available for the workers.
            assignin('base', 'Fs', Fs);
            assignin('base', 'USRP_NOISE_FLOOR_V', USRP_NOISE_FLOOR_V);
            curRx3D = curPathLossUtmXYHs(idxSeg, :);
            curRx3D(:,3) = curRx3D(:,3) + RX_HEIGHT_M;
            
            curDist = norm([xTx, yTx, TX_ALT+TX_HEIGHT_M] - curRx3D);
            
            if(curDist>=curRange(1) && curDist<=curRange(2)) %#ok<PFBNS>
                % Plot and save the PDP.
                curSigOutFileRelPath ...
                    = curContiOutFilesRelPathsUnderDataFolder;
                curSigOutFile = rdir(...
                    fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'Data', ...
                    curSigOutFileRelPath));
                
                curSegRange = contiPathLossesExtraInfo.contiSampIndices...
                    {idxTrack}(idxSeg,:);
                
                [hPdpsFig, timeMsInPlot, signalAmpsInPlot] ...
                    = plotPdpsForOneRec(curSigOutFile, F_S, curSegRange);
                
                plotFileName = [...
                    'PdpOverview_Dist_', num2str(curDist, '%.0f'), ...
                    '_Track_', num2str(idxTrack), ...
                    '_Seg_', num2str(idxSeg)];
                saveas(hPdpsFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
                    [plotFileName, '.jpg']));
                
                % Also estimate the energy ratio of the LoS peak.
                losPeakEnergyRsSegs(idxSeg) ...
                    = estimateEnergyRatioInOnePdpForLosSig( ...
                    timeMsInPlot, signalAmpsInPlot, ...
                    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
                    [plotFileName, '_peaks.jpg']));
                
                % Close the figure.
                close(hPdpsFig);
            end
        end
        
        losPeakEnergyRsTracks{idxTrack} = losPeakEnergyRsSegs;
    end
    
    losPeakEnergyRs{idxRange} = losPeakEnergyRsTracks;
end

% Save and plot the distribution of the LoS peak energy ratios.
save(fullfile(ABS_PATH_TO_SAVE_PLOTS, 'losPeakEnergyRs.mat'), ...
    'losPeakEnergyRs');
for idxRange = 1:numDistRangesToInsp
    curRange = distRangesToInsp{idxRange};
    % Get all the energy ratio results (from all segments of different
    % tracks) within that range.
    curPeakEngergyRs = vertcat(losPeakEnergyRs{idxRange}{:});
    
    % Empirical CDF.
    [r,x] = ecdf(curPeakEngergyRs);
    
    hPeakEnRsFig = figure; hold on;    
    plot(x, r, 'b.-'); 
    xlabel('LoS Peak Energy Ratio'); ylabel('Empirical CDF');
    plotFileName = [...
        'losPeakEnergyRs_Range_', num2str(curRange(1)), ...
        '_to_', num2str(curRange(2))];    
    saveas(hPeakEnRsFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [plotFileName, '.jpg']));
end

disp('    Done!')

% EOF