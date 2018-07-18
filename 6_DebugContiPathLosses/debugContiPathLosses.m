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

% For loading TX location.
addpath(fullfile(pwd, '2_PlotGpsSampsOnMap'));
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

% Generate debug figures for noise elimination or not.
flagGenerateNoiseEliDebugFig = true;

% For generating .eps figures.
pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', '1_EpsFigs');

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

% Delete all files in under ABS_PATH_TO_SAVE_PLOTS.
dinfo = dir(ABS_PATH_TO_SAVE_PLOTS);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(ABS_PATH_TO_SAVE_PLOTS, {dinfo.name});
warning('We need to delete history results from this script!')
k = input('    Is this OK? (y/n) ', 's');
if ~strcmpi(k, 'y')
    error('User declined to continue.')
end
if ~isempty(filenames)
    delete( filenames{:} )
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
% NUM_SIGMA_FOR_THRESHOLD = 3;

disp('    Done!')

%% Inspect PDP's within Specified Diantance Ranges

disp('    Generating PDP plots for signal recording segments ...');

% Set this to control how to group the sample locations according to the
% RX-to-TX distance.
%   For example, {[50, 500], [0, 50], [175, 180]}: The first range contains
%	almost everything; Other ranges are for separate analyses on smaller
%	ranges.
distRangesToInsp = {[0, 500]}; % In meter.
numDistRangesToInsp = length(distRangesToInsp);

losPeakEnergyRs = cell(numDistRangesToInsp, 1);
% We need (lat, lon) for the segments for illustrating results on a map.
losPeakEstiLatLon = cell(numDistRangesToInsp, 1);
for idxRange = 1:numDistRangesToInsp
    curRange = distRangesToInsp{idxRange};
    curNumTracks = length(contiPathLossesWithGpsInfo);
    
    losPeakEnergyRsTracks = cell(curNumTracks, 1);
    estiLatLonSegsTracks = cell(curNumTracks, 1);
    
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
        curContiPathLossesWithGpsInfo ...
            = contiPathLossesWithGpsInfo{idxTrack};
        
        parfor idxSeg = 1:numCurSegs
            % Make necessary variables available for the workers.
            assignin('base', 'Fs', Fs); %#ok<PFEVB>
            assignin('base', ...
                'USRP_NOISE_FLOOR_V', USRP_NOISE_FLOOR_V); %#ok<PFEVB>
            assignin('base', 'FLAG_PDP_TIME_REVERSED', ...
                FLAG_PDP_TIME_REVERSED); %#ok<PFEVB>
            assignin('base', 'SLIDE_FACTOR', SLIDE_FACTOR); %#ok<PFEVB>
            assignin('base', ...
                'NUM_SIGMA_FOR_THRESHOLD', ...
                NUM_SIGMA_FOR_THRESHOLD); %#ok<PFEVB>
            
            curRx3D = curPathLossUtmXYHs(idxSeg, :);
            curRx3D(:,3) = curRx3D(:,3) + RX_HEIGHT_M;
            
            curDist = norm([xTx, yTx, TX_ALT+TX_HEIGHT_M] - curRx3D);
            
            if(curDist>=curRange(1) && curDist<=curRange(2))
                % Plot and save the PDP.
                curSigOutFileRelPath ...
                    = curContiOutFilesRelPathsUnderDataFolder;
                curSigOutFile = rdir(...
                    fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'Data', ...
                    curSigOutFileRelPath));
                
                curSegRange = contiPathLossesExtraInfo.contiSampIndices...
                    {idxTrack}(idxSeg,:);
                
                [hPdpsFig, timeMsInPlot, signalAmpsInPlot, ...
                    lowPassedSigInPlot, ~, hNoiseEliDebugFig] ...
                    = plotPdpsForOneRec(curSigOutFile, F_S, ...
                    curSegRange, flagGenerateNoiseEliDebugFig);
                
                plotFileName = [...
                    'PdpOverview_Dist_', num2str(curDist, '%.0f'), ...
                    '_Track_', num2str(idxTrack), ...
                    '_Seg_', num2str(idxSeg)];
                saveas(hPdpsFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
                    [plotFileName, '.jpg']));
                % Close the figure.
                close(hPdpsFig);
                
                % Also estimate the energy ratio of the LoS peak.
                losPeakEnergyRsSegs(idxSeg) ...
                    = estimateEnergyRatioInOnePdpForLosSig( ...
                    timeMsInPlot, signalAmpsInPlot, ...
                    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
                    [plotFileName, '_Peaks.jpg']), lowPassedSigInPlot);
                
                % Save the noise elimination figure.
                if isgraphics(hNoiseEliDebugFig) ...
                        && isvalid(hNoiseEliDebugFig)
                    saveas(hNoiseEliDebugFig, ...
                        fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
                        [plotFileName, '_DynamicNoiseEli.jpg']));
                    close(hNoiseEliDebugFig);
                end
            end
        end
        
        losPeakEnergyRsTracks{idxTrack} = losPeakEnergyRsSegs;
        estiLatLonSegsTracks{idxTrack} ...
            = curContiPathLossesWithGpsInfo(:, 2:3);
    end
    
    losPeakEnergyRs{idxRange} = losPeakEnergyRsTracks;
    losPeakEstiLatLon{idxRange} = estiLatLonSegsTracks;
end

% Save and plot the distribution of the LoS peak energy ratios.
save(fullfile(ABS_PATH_TO_SAVE_PLOTS, 'losPeakEnergyRs.mat'), ...
    'losPeakEnergyRs', 'losPeakEstiLatLon');

%% Generate an Overview Figure
% Load the TX location.
[markerLats, markerLons, markerNames] ...
    = loadGpsMarkers(absPathToGpsMarkerCsvFile);
idxTxMarker = find(strcmp(markerNames, 'Tx'));
latTx = markerLats(idxTxMarker);
lonTx = markerLons(idxTxMarker);
for idxRange = 1:numDistRangesToInsp
    curRange = distRangesToInsp{idxRange};
    % Get all the energy ratio results (from all segments of different
    % tracks) within that range.
    curLosPeakEngergyRs = vertcat(losPeakEnergyRs{idxRange}{:});
    
    % Empirical CDF.
    [r,x] = ecdf(curLosPeakEngergyRs);
    
    hPeakEnRsCdfFig = figure; hold on;
    plot(x, r, '--', 'Color', ones(1,3)*0.5);
    plot(x, r, 'b.');
    xlabel('LoS Peak Energy Ratio'); ylabel('Empirical CDF'); grid on;
    plotFileName = [...
        'losPeakEnergyRs_CDF_Range_', num2str(curRange(1)), ...
        '_to_', num2str(curRange(2))];
    saveas(hPeakEnRsCdfFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [plotFileName, '.fig']));
    saveas(hPeakEnRsCdfFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [plotFileName, '.jpg']));
    
    % Export an .eps copy for papers.
    saveEpsFigForPaper(hPeakEnRsCdfFig, ...
        fullfile(pathToSavePaperFigs, ...
        ['5_', plotFileName, '.eps']));
    
    % Plot the LoS peak energy ratio results on a map to show their
    % locations.
    curLosPeakEstiLatLon = vertcat(losPeakEstiLatLon{idxRange}{:});
    
    hPeakEnRsOnMapFig = figure; hold on;
    boolsNanLosPeakRs = isnan(curLosPeakEngergyRs);
    % Valid results.
    %     colormap('jet');
    plot3k([curLosPeakEstiLatLon(~boolsNanLosPeakRs,2), ...
        curLosPeakEstiLatLon(~boolsNanLosPeakRs,1), ...
        curLosPeakEngergyRs(~boolsNanLosPeakRs)]);
    plot3(curLosPeakEstiLatLon(boolsNanLosPeakRs,2), ...
        curLosPeakEstiLatLon(boolsNanLosPeakRs,1), ...
        curLosPeakEngergyRs(boolsNanLosPeakRs), 'xk');
    hTx = plot3(lonTx, latTx, 0, 'g^');
    legend(hTx, 'Tx'); xlabel('Longitude'); ylabel('Latitude');
    xticks([]); yticks([]); grid on; view(2);
    
    % Manually adjust the visible area.
    axis([-105.2776 -105.2743 39.9892 39.9917]);
    plot_google_map('MapType', 'satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = findall( allchild(hPeakEnRsOnMapFig), 'type', 'colorbar');
    hCbTicksToSet = 0:0.1:1;
    hCb.TickLabels = arrayfun(@(n) {num2str(n)}, hCbTicksToSet)';
    hCb.Ticks = linspace(1,length(colormap)+1, ...
        length(hCb.TickLabels));
    
    plotFileName = [...
        'losPeakEnergyRs_map_Range_', num2str(curRange(1)), ...
        '_to_', num2str(curRange(2))];
    saveas(hPeakEnRsOnMapFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [plotFileName, '.fig']));
    saveas(hPeakEnRsOnMapFig, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [plotFileName, '.jpg']));
end

disp('    Done!')

% EOF