% EXPORTDATAFORSIMULATION Export data as csv files for external usages,
% e.g. simulations.
%
% Yaguang Zhang, Purdue, 12/27/2018

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% We will need the functions latLon2PixIndices.m and pixIndices2LatLon.m
% for working with the vegArea image.
addpath(fullfile(pwd, '9_GenerateVegAreas'));

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'SimulationDataExportation');

% Reuse results from evalPathLossesForContiTracks.m and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

% Reuse results from generateVegAreas.m.
ABS_PATH_TO_VEG_AREAS_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'AutoGeneratedVegAreas', 'vegAreasMeta.mat');

% Reuse results from estimateFoliageAttenuation.m.
ABS_PATH_TO_UTM_INFO_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'FoliageAttenuationEstimation_ManualTreeLocs', ...
    'utmInfoForPathLossesAndTrees.mat');

% Reuse results from estimateFoliageAttenuationWithManualTreeLocs.m.
ABS_PATH_TO_TREE_NUM_BASED_ANALYSIS_FILE ...
    = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'PostProcessingResults', ...
    'FoliageAttenuationEstimation_ManualTreeLocs', ...
    'foliageAttenAnalysisResults.mat');

% Control the truncation of the data in the output .csv files.
dataStrFormatter = '%.12f';

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
disp('      - vegAreasMeta.mat')
disp('      - utmInfoForPathLossesAndTrees.mat')

assert(exist(ABS_PATH_TO_PATH_LOSSES_FILE, 'file')==2, ...
    'Couldn''t find contiPathLossesWithGpsInfo.mat! Please run NistMeasurementCampaignCode/3_PathLossComputation/evalPathLossesForContiTracks.m first.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run NistMeasurementCampaignCode/3_PathLossComputation/loadMeasCampaignInfo.m first.');
assert(exist(ABS_PATH_TO_VEG_AREAS_FILE, 'file')==2, ...
    'Couldn''t find vegAreasMeta.mat! Please run NistMeasurementCampaignCode/9_GenerateVegAreas/generateVegAreas.m first.');
assert(exist(ABS_PATH_TO_UTM_INFO_FILE, 'file')==2, ...
    'Couldn''t find utmInfoForPathLossesAndTrees.mat! Please run NistMeasurementCampaignCode/4_FoliageAttenuationEstimation/estimateFoliageAttenuation.m first.');
assert(exist(ABS_PATH_TO_TREE_NUM_BASED_ANALYSIS_FILE, 'file')==2, ...
    'Couldn''t find foliageAttenAnalysisResults.mat! Please run NistMeasurementCampaignCode/8_FoliageAttenuationEstimation_ManualTreeLocs/estimateFoliageAttenuationWithManualTreeLocs.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'contiPathLossesWithGpsInfo', 'contiOutFilesRelPathsUnderDataFolder'
% and 'contiOutFileIndicesReflection'.
load(ABS_PATH_TO_PATH_LOSSES_FILE);
% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 	'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Get 'vegAreas', 'LAT_RANGE', 'LON_RANGE', 'VEG_AREA_IMG_RESOLUTION', and
% 'VEG_AREA_IMG_META'.
load(ABS_PATH_TO_VEG_AREAS_FILE);
% Get 'xTx', 'yTx', 'txUtmZone', 'treeLocations', 'pathLossUtmXYHs',
% 'pathLossUtmZones', 'treeUtmXYHs', and 'treeUtmZones'.
load(ABS_PATH_TO_UTM_INFO_FILE);
% Get 'numsOfTreesInFirstFresnel', 'freeSpacePathLosses',
% 'exceLossRefFreeSpace', 'excePathLossPerTree',
% 'shiftedFreeSpacePathLosses', 'excePathLossGroupWise', and
% 'groupWiseShiftedFreeSpacePathLosses'.
load(ABS_PATH_TO_TREE_NUM_BASED_ANALYSIS_FILE);

disp('    Done!')

%% Export Trunk Locations

% Shift negative tree height values for the trunk locations to this.
MIN_TREE_HEIGHT_ALLOWED_IN_METER = 0.25;

disp('Collecting foliage area data ...')
groundHeightWrtTXInM = VEG_AREA_IMG_META.ALTS - (TX_ALT+TX_HEIGHT_M);
treeHeightInM = VEG_AREA_IMG_META.ZS - VEG_AREA_IMG_META.ALTS;

getInterGroundHeightWrtTXInM = scatteredInterpolant( ...
    VEG_AREA_IMG_META.XS(:), ...
    VEG_AREA_IMG_META.YS(:), ...
    groundHeightWrtTXInM(:));
getInterTreeHeightInM = scatteredInterpolant( ...
    VEG_AREA_IMG_META.XS(:), ...
    VEG_AREA_IMG_META.YS(:), ...
    treeHeightInM(:));

trunkGroundHeightWrtTXInM ...
    = getInterGroundHeightWrtTXInM(treeUtmXYHs(:,1), treeUtmXYHs(:,2));
trunkTreeHeightInM ...
    = getInterTreeHeightInM(treeUtmXYHs(:,1), treeUtmXYHs(:,2));

% Fix the negative tree height issue.
boolsTrunkLocsWithNegTreeH = trunkTreeHeightInM<0;
trunkTreeHeightInM(boolsTrunkLocsWithNegTreeH) ...
    = MIN_TREE_HEIGHT_ALLOWED_IN_METER;

% Ignore entries with NaN attribute(s).
boolsValidEntries ...
    = (~isnan(trunkGroundHeightWrtTXInM))...
    &(~isnan(trunkTreeHeightInM));

fullPathTrunkLocCsv = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocs.csv');
curHeaderCell = {'utmX', 'utmY', 'utmZone', ...
    'lat', 'lon', 'groundHeightWrtTXInM', 'treeHeightInM'};
curData = [num2cell(treeUtmXYHs(boolsValidEntries,1:2)), ...
    cellstr(treeUtmZones(boolsValidEntries,:)), ...
    num2cell(treeLocations(boolsValidEntries,1:2)), ...
    num2cell(trunkGroundHeightWrtTXInM(boolsValidEntries, :)), ...
    num2cell(trunkTreeHeightInM(boolsValidEntries, :))];
disp('Done!')

disp('Writing trunk locations to file ...')
writeToCsvWithHeader(fullPathTrunkLocCsv, curHeaderCell, curData, dataStrFormatter);
disp('Done!')

% Illustration plots.
curFigTitle = 'Tree Locations - (lon, lat)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsLonLat.png');

curFig = figure;
plot(treeLocations(:,2), treeLocations(:,1), 'rx');
plot_google_map('MapType', 'satellite');
xlabel('lon'); ylabel('lat');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = 'Tree Locations - (utmX, utmY)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsUtmXY.png');

curFig = figure;
plot(treeUtmXYHs(:,1), treeUtmXYHs(:,2), 'rx');
grid on; axis equal;
xlabel('utmX'); ylabel('utmY');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = 'Tree Locations with Tree Height - (utmX, utmY, treeHeightInM)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsUtmXYTreeHeight.png');

curFig = figure;
plot3k([treeUtmXYHs(:,1), treeUtmXYHs(:,2), trunkTreeHeightInM]);
grid on; axis equal;
xlabel('utmX'); ylabel('utmY');

title(curFigTitle)
saveas(curFig, curFigPath);

%% Export Foliage Areas

fullPathFoliageAreaCsv = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'foliageAreaGrid.csv');
curHeaderCell = {'rowIndexM', 'columnIndexN', ...
    'utmX', 'utmY', 'utmZone', 'lat', 'lon', 'altInMFromNist'...
    'groundHeightWrtTXInM', 'treeHeightInM', ...
    'booleanInFoliage'};

numRows = VEG_AREA_IMG_META.IMG_RESOLUTION(2);
numCols = VEG_AREA_IMG_META.IMG_RESOLUTION(1);

% Because of the huge data size, we need to write into the file in batches,
% instead of:
%   writeToCsvWithHeader(fullPathFoliageAreaCsv, curHeaderCell, curData);
curData = cell(numCols, length(curHeaderCell));

% File header.
cell2csv(fullPathFoliageAreaCsv, curHeaderCell);

disp('Collecting foliage area data ...')
for m = 1:numRows
    disp(['    Row ', num2str(m), '/' num2str(numRows)]);
    
    dataRowCnt = 1;
    for n = 1:numCols
        curData(dataRowCnt,:) = {m, n, ...
            VEG_AREA_IMG_META.XS(m,n), VEG_AREA_IMG_META.YS(m,n), ...
            VEG_AREA_IMG_META.UTM_ZONE, ...
            VEG_AREA_IMG_META.LATS(m,n), VEG_AREA_IMG_META.LONS(m,n), ...
            VEG_AREA_IMG_META.ALTS(m,n), ...
            groundHeightWrtTXInM(m,n), treeHeightInM(m,n), ...
            vegAreas(m,n)};
        dataRowCnt = dataRowCnt+1;
    end
    appendCell2Csv(fullPathFoliageAreaCsv, curData, dataStrFormatter);
end
disp('Done!')

% Illustration plots.
curFigTitle = {'Ground Height Wrt TX', 'for the Foliage Areas'};
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'foliageAreasGroundHeightWrtTX.png');

curFig = figure;
plot3k([VEG_AREA_IMG_META.XS(vegAreas(:)==1), ...
    VEG_AREA_IMG_META.YS(vegAreas(:)==1), ...
    groundHeightWrtTXInM(vegAreas(:)==1)]);
grid on; axis equal; view(2);
xlabel('utmX'); ylabel('utmY'); zlabel('groundHeightWrtTXInM');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = {'Tree Height for the Foliage Areas'};
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'foliageAreasTreeHeight.png');

curFig = figure;
plot3k([VEG_AREA_IMG_META.XS(vegAreas(:)==1), ...
    VEG_AREA_IMG_META.YS(vegAreas(:)==1), ...
    treeHeightInM(vegAreas(:)==1)]);
grid on; axis equal; view(2);
xlabel('utmX'); ylabel('utmY'); zlabel('treeHeightInM');

title(curFigTitle)
saveas(curFig, curFigPath);

%% Export TX and RX Information

fullPathAntennaAngles = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'antennaAngles.csv');

numTracks = length(TX_INFO_LOGS{1});
curData = struct();
curData.trackIndex = (1:numTracks)';
curData.txAzimuth = {TX_INFO_LOGS{1}.txAz}';
curData.txElevation = {TX_INFO_LOGS{1}.txEl}';
curData.rxAzimuth = {TX_INFO_LOGS{1}.rxAz}';
curData.rxElevation = {TX_INFO_LOGS{1}.rxEl}';
curData.rxHeightAboveGroundInM = 1;
curData.description = {TX_INFO_LOGS{1}.location}';

% Save to file.
struct2csv(curData, fullPathAntennaAngles);
% Remove double quotes and tailing comma.
curDataStr = fileread(fullPathAntennaAngles);
curDataStr = erase(curDataStr, '"');
curDataStr = regexprep(curDataStr, ',\n', '\n');

curFileHandle = fopen(fullPathAntennaAngles, 'w');
fwrite(curFileHandle, curDataStr, 'char');
fclose(curFileHandle);

fullPathTxLoc = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'txLoc.csv');
curHeaderCell = {'utmX', 'utmY', 'utmZone', ...
    'lat', 'lon', 'altFromGoogle', ...
    'groundHeightWrtTXInM', ...
    'heightAboveGroundInFeet', 'heightAboveGroundInM'};
curData = {xTx, yTx, txUtmZone, ...
    TX_LAT, TX_LON, TX_ALT, ...
    -TX_HEIGHT_M, ...
    TX_HEIGHT_FEET, TX_HEIGHT_M};
writeToCsvWithHeader(fullPathTxLoc, curHeaderCell, curData, dataStrFormatter);

txInfoLog = TX_INFO_LOGS{1};
for idxTrack = 1:numTracks
    % RX location.
    fullPathRxLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['rxLoc_track_', num2str(idxTrack), '.csv']);
    curHeaderCell = {'utmX', 'utmY', 'utmZone', ...
        'lat', 'lon', 'altFromGoogle', 'rxHeightWrtTXInM'};
    curXsRx = pathLossUtmXYHs{idxTrack}(:,1);
    curYsRx = pathLossUtmXYHs{idxTrack}(:,2);
    curRxHeightWrtTxInM = (RX_HEIGHT_M ...
        +contiPathLossesWithGpsInfo{idxTrack}(:,4))...
        -(TX_ALT+TX_HEIGHT_M);
    curData = [num2cell([curXsRx, ...
        curYsRx]), ...
        cellstr(pathLossUtmZones{idxTrack}), ...
        ...
        num2cell([contiPathLossesWithGpsInfo{idxTrack}(:,2), ...
        contiPathLossesWithGpsInfo{idxTrack}(:,3), ...
        contiPathLossesWithGpsInfo{idxTrack}(:,4), ...
        curRxHeightWrtTxInM])];
    writeToCsvWithHeader(fullPathRxLocs, curHeaderCell, curData, dataStrFormatter);
    
    % Measurement results.
    fullPathRxMeas = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['rxLoc_meas_', num2str(idxTrack), '.csv']);
    curHeaderCell = {'locIdx', 'pathLossInDb', 'rxToTx3DDistInM', ...
        'txAzimuth (clockwise from positive y)', ...
        'txElevation ("+" for upward)', ...
        'rxAzimuth (clockwise from positive y)', ...
        'rxElevation ("+" for upward)'};
    curMeas = contiPathLossesWithGpsInfo{idxTrack}(:,1);
    curRxToTx3DDistInM = vecnorm([curXsRx - xTx, ...
        curYsRx - yTx, curRxHeightWrtTxInM], 2, 2);
    
    onesWithNumOfCurMeas = ones(length(curMeas), 1);
    curData = [(1:length(curMeas))', curMeas, curRxToTx3DDistInM, ...
        onesWithNumOfCurMeas.*(360 - txInfoLog(idxTrack).txAz), ...
        onesWithNumOfCurMeas.*txInfoLog(idxTrack).txEl, ...
        onesWithNumOfCurMeas.*(360 - txInfoLog(idxTrack).rxAz), ...
        onesWithNumOfCurMeas.*txInfoLog(idxTrack).rxEl];
    writeToCsvWithHeader(fullPathRxMeas, curHeaderCell, curData, dataStrFormatter);
end

%% Illustration plots.
curFigTitle = 'TX and RX Locations - (lon, lat)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'txAndRxLocsLonLat.png');

curFig = figure; hold on;
hTx = plot(TX_LON, TX_LAT, 'gx');
for idxTrack = 1:numTracks
    hRx = plot(contiPathLossesWithGpsInfo{idxTrack}(:,3), ...
        contiPathLossesWithGpsInfo{idxTrack}(:,2), 'r.', 'MarkerSize', 5);
end
axis tight;
plot_google_map('MapType', 'satellite');
xlabel('lon'); ylabel('lat');
legend([hTx, hRx], 'TX', 'RX', 'Location', 'SouthEast');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = 'TX and RX Locations - (utmX, utmY)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'txAndRxLocsUtmXY.png');

curFig = figure; hold on;
hTx = plot(xTx, yTx, 'gx');
for idxTrack = 1:numTracks
    hRx = plot(pathLossUtmXYHs{idxTrack}(:,1), ...
        pathLossUtmXYHs{idxTrack}(:,2), 'r.', 'MarkerSize', 5);
end
grid on; axis equal;
xlabel('utmX'); ylabel('utmY');
legend([hTx, hRx], 'TX', 'RX', 'Location', 'SouthEast');

title(curFigTitle)
saveas(curFig, curFigPath);

%% Extra Manually Generated Plots (Not Saved Automatically)
FLAG_GENERATE_EXTRA_ILLUSTRATIONS = false;

if FLAG_GENERATE_EXTRA_ILLUSTRATIONS
    % Illustration for the Grid Points Colored by LiDAR Z Values
    xs = VEG_AREA_IMG_META.XS(1:300, 1:300);
    ys = VEG_AREA_IMG_META.YS(1:300, 1:300);
    zs = VEG_AREA_IMG_META.ZS(1:300, 1:300);
    figure('position', [0 0 700 500]); plot3k([xs(:), ys(:), zs(:)])
    axis([min(xs(:)) max(xs(:)) min(ys(:)) max(ys(:)) min(zs(:)) max(zs(:))])
    xlabel('UTM x (m)'); ylabel('UTM y (m)');
    title({'Illustration for the Grid Points';'Colored by LiDAR Z Values'})
    view(2)
    
    % Foliage Area Colored by Tree Height with TX, RX and Trunk Locations
    figure('position', [0 0 700 500]); hold on;
    foliageAreaTreeHeightInM = treeHeightInM(vegAreas(:)==1);
    plot3k([VEG_AREA_IMG_META.XS(vegAreas(:)==1), ...
        VEG_AREA_IMG_META.YS(vegAreas(:)==1), ...
        foliageAreaTreeHeightInM]);
    axis([min(VEG_AREA_IMG_META.XS(:)) max(VEG_AREA_IMG_META.XS(:)) ...
        min(VEG_AREA_IMG_META.YS(:)) max(VEG_AREA_IMG_META.YS(:)) ...
        min(foliageAreaTreeHeightInM) max(foliageAreaTreeHeightInM)])
    xlabel('UTM x (m)'); ylabel('UTM y (m)');
    hTx = plot3(xTx, yTx, max(foliageAreaTreeHeightInM), 'go', ...
        'LineWidth', 2);
    for idxTrack = 1:length(TX_INFO_LOGS{1})
        hRx = plot3(pathLossUtmXYHs{idxTrack}(:,1), ...
            pathLossUtmXYHs{idxTrack}(:,2), ...
            ones(length(pathLossUtmXYHs{idxTrack}(:,1)),1)...
            .*max(foliageAreaTreeHeightInM), ...
            'r.', 'MarkerSize', 5);
    end
    hTrunks = plot3(treeUtmXYHs(:,1), treeUtmXYHs(:,2), ...
        ones(length(treeUtmXYHs(:,1)),1) ...
        .*max(foliageAreaTreeHeightInM), 'kx', 'LineWidth', 2);
    legend([hTx, hRx, hTrunks], 'TX', 'RX', 'Trunks', ...
        'Location', 'SouthEast');
    title({'Foliage Area Colored by Tree Height'; ...
        'with TX, RX and Trunk Locations'})
    view(2)
    
    % Data availability inspection.
    BUFFER_DIST_IN_M = 25;
    alphaForTranparentObjs = 0.3;
    figure('position', [0 0 700 500]); hold on;
    xRec = min(VEG_AREA_IMG_META.XS(:));
    yRec = min(VEG_AREA_IMG_META.YS(:));
    wRec = max(VEG_AREA_IMG_META.XS(:)) - min(VEG_AREA_IMG_META.XS(:));
    hRec = max(VEG_AREA_IMG_META.YS(:)) - min(VEG_AREA_IMG_META.YS(:));
    hFoliageAreaData = rectangle('Position', [xRec, yRec, wRec, hRec], ...
        'FaceColor', [ones(1,3).*0.8 alphaForTranparentObjs]);
    for idxTrack = 1:length(TX_INFO_LOGS{1})
        numRxLocs = length(pathLossUtmXYHs{idxTrack}(:,1));
        hAreaNeeded = viscircles( ...
            [pathLossUtmXYHs{idxTrack}(:,1) ...
            pathLossUtmXYHs{idxTrack}(:,2)], ...
            ones(numRxLocs,1).*BUFFER_DIST_IN_M, 'color', 'y');
    end
    hTx = plot(xTx, yTx, 'go', ...
        'LineWidth', 2);
    for idxTrack = 1:length(TX_INFO_LOGS{1})
        hRx = plot(pathLossUtmXYHs{idxTrack}(:,1), ...
            pathLossUtmXYHs{idxTrack}(:,2), ...
            'r.', 'MarkerSize', 5);
    end
    hTrunks = plot(treeUtmXYHs(:,1), treeUtmXYHs(:,2),...
        'kx', 'LineWidth', 2);
    axis equal;
    xlabel('UTM x (m)'); ylabel('UTM y (m)');
    legend([hTx, hRx, hTrunks, hAreaNeeded], ...
        'TX', 'RX', 'Trunks', ...
        [num2str(BUFFER_DIST_IN_M), ' m buffer area'],...
        'Location', 'SouthEast');
    title('Data availability inspectation')
    
    % Inspect the negtive tree height issue.
    indicesForTrunksWithNegHeight = trunkTreeHeightInM<0;
    figure('position', [0 0 700 500]);
    plot(treeLocations(indicesForTrunksWithNegHeight,2), ...
        treeLocations(indicesForTrunksWithNegHeight,1), 'rx');
    plot_google_map('MapType', 'satellite');
    xticks([]); yticks([]);
    title('Trunk Locations with Negtive Tree Height');
end

%% Close Figures
close all;

% EOF