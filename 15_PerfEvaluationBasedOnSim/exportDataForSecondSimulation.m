% EXPORTDATAFORSIMULATION Export data as csv files for external usages,
% e.g. simulations.
%
% We will export:
%   1. TX information;
%    2. Measurement RX location informaiton;
%   3. Relative geographical data, including elevation and LiDAR
%   information;
%    4. Tree locations and their heights.
%
% Yaguang Zhang, Purdue, 09/23/2019

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
    'PostProcessingResults', 'SimulationDataExportationExtended');

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

% The updated tree locations.
ABS_PATH_TO_MANUALLY_LABELED_TREE_LOCS ...
    = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'PostProcessingResults', ...
    'ManuallyLocateMoreTrees', 'treeLocs.mat');

% Control the truncation of the data in the output .csv files.
dataStrFormatter = '%.12f';

% For plotting.
lightGray = ones(1,3).*0.7;

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
disp('      - treeLocs.mat (Extended)')

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
assert(exist(ABS_PATH_TO_MANUALLY_LABELED_TREE_LOCS, 'file')==2, ...
    'Couldn''t find treeLocs.mat! Please run NistMeasurementCampaignCode/15_PerfEvaluationBasedOnSim/manuallyLocateMoreTrees.m first.');

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
% Get 'markLocs'.
load(ABS_PATH_TO_MANUALLY_LABELED_TREE_LOCS);

disp('    Done!')

%% Interpolate Data in the Vegetation Grid

% Shift negative tree height values for the trunk locations to this.
MIN_TREE_HEIGHT_ALLOWED_IN_METER = 0.25;

disp(' ')
disp('    Collecting foliage area data ...')

% Function to fetch elevation according to the vegetation data structure.
getUsgsAltInM = scatteredInterpolant( ...
    VEG_AREA_IMG_META.XS(:), ...
    VEG_AREA_IMG_META.YS(:), ...
    VEG_AREA_IMG_META.ALTS(:));
TX_ALT_USGS = getUsgsAltInM(xTx, yTx);

groundHeightWrtTXInM = VEG_AREA_IMG_META.ALTS - (TX_ALT_USGS+TX_HEIGHT_M);
treeHeightInM = VEG_AREA_IMG_META.ZS - VEG_AREA_IMG_META.ALTS;

% Grid data for ground height w.r.t. the TX antenna.
getInterGroundHeightWrtTXInM = @(xs, ys) ...
    getUsgsAltInM(xs,ys) - (TX_ALT_USGS+TX_HEIGHT_M);
% Grid data for vegetation height.
getLidarZInM = scatteredInterpolant( ...
    VEG_AREA_IMG_META.XS(:), ...
    VEG_AREA_IMG_META.YS(:), ...
    VEG_AREA_IMG_META.ZS(:));
getInterTreeHeightInM = @(xs, ys) ...
    getLidarZInM(xs,ys) - getUsgsAltInM(xs,ys);

disp('Done!')

%%  Export Trunk Locations

% Note: "treeUtmXYsExt" contains the previous trunk location logs for the
% ICC wireless letter paper; here we will use an extended set of trunk
% location logs, named "treeUtmXYsExt", which contains more trees to cover
% the whole LiDAR grid area.

disp(' ')
disp('    Fetching tree locations ...')

% Remove trees that are out of the LiDAR coverage area.
minLidarLat = min(VEG_AREA_IMG_META.LATS(:));
maxLidarLat = max(VEG_AREA_IMG_META.LATS(:));
minLidarLon = min(VEG_AREA_IMG_META.LONS(:));
maxLidarLon = max(VEG_AREA_IMG_META.LONS(:));

boolsOutOfLidarCovArea = markLocs(:,1)>=minLidarLat ...
    & markLocs(:,1)<=maxLidarLat ...
    & markLocs(:,2)>=minLidarLon ...
    & markLocs(:,2)<=maxLidarLon;
if ~exist('markLocsOriginal', 'var')
    markLocsOriginal = markLocs;
    [numTreesOriginal, ~] = size(markLocs);
end
markLocs(~boolsOutOfLidarCovArea,:) = [];

% Covert the GPS tree locations to UTM.
[numTrees, ~] = size(markLocs); treeUtmXYsExt = nan(numTrees, 2);
[treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), treeUtmZonesExt] ...
    = deg2utm(markLocs(:,1), markLocs(:,2));

assert(all([arrayfun(@(idx) ...
    strcmp(treeUtmZonesExt(idx,:), treeUtmZonesExt(1,:)), ...
    1:numTrees)]), 'All trees should be in the same UTM zone!');

trunkGroundHeightWrtTXInM ...
    = getInterGroundHeightWrtTXInM(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2));
trunkTreeHeightInM ...
    = getInterTreeHeightInM(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2));

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
curData = [num2cell(treeUtmXYsExt(boolsValidEntries,1:2)), ...
    cellstr(treeUtmZonesExt(boolsValidEntries,:)), ...
    num2cell(markLocs(boolsValidEntries,1:2)), ...
    num2cell(trunkGroundHeightWrtTXInM(boolsValidEntries, :)), ...
    num2cell(trunkTreeHeightInM(boolsValidEntries, :))];

disp('Done!')

disp(' ')
disp('Writing trunk locations to file ...')
writeToCsvWithHeader(fullPathTrunkLocCsv, curHeaderCell, ...
    curData, dataStrFormatter);
disp('Done!')

% Illustration plots.
curFigTitle = 'Tree Locations - (lon, lat)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsLonLat.png');

curFig = figure; hold on;
hOriTrees = plot(markLocsOriginal(:,2), markLocsOriginal(:,1), 'r.');
hTrees = plot(markLocs(boolsValidEntries,2), ...
    markLocs(boolsValidEntries,1), 'b*');
plot_google_map('MapType', 'satellite');
plot3(markLocsOriginal(:,2), markLocsOriginal(:,1), ...
    ones(numTreesOriginal,1).*(max(VEG_AREA_IMG_META.ZS)+1),  'r.');
plot3(markLocs(boolsValidEntries,2), ...
    markLocs(boolsValidEntries,1), ...
    ones(sum(boolsValidEntries),1).*(max(VEG_AREA_IMG_META.ZS)+1), 'b*');
plot3k([VEG_AREA_IMG_META.LONS(:), VEG_AREA_IMG_META.LATS(:), ...
    VEG_AREA_IMG_META.ZS(:)], 'ColorBar', false);
xlabel('lon'); ylabel('lat'); view(2);
legend([hOriTrees, hTrees], 'All trees', 'Valid trees', ...
    'Location', 'southeast');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = 'Tree Locations - (utmX, utmY)';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsUtmXY.png');

curFig = figure;
plot(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), 'rx');
grid on; axis equal;
xlabel('utmX'); ylabel('utmY');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle ...
    = 'Tree Locations with Tree Height - (utmX, utmY, treeHeightInM)';
curFigPath ...
    = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'trunkLocsUtmXYTreeHeight.png');

curFig = figure;
plot3k([treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), trunkTreeHeightInM]);
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

disp(' ')
disp('Exporting foliage area data ...')
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

fullPathAntennaAngles ...
    = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'antennaAngles.csv');

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
    'lat', 'lon', 'altFromUsgs', ...
    'groundHeightWrtTXInM', ...
    'heightAboveGroundInFeet', 'heightAboveGroundInM'};
curData = {xTx, yTx, txUtmZone, ...
    TX_LAT, TX_LON, TX_ALT_USGS, ...
    -TX_HEIGHT_M, ...
    TX_HEIGHT_FEET, TX_HEIGHT_M};
writeToCsvWithHeader(fullPathTxLoc, curHeaderCell, curData, dataStrFormatter);

txInfoLog = TX_INFO_LOGS{1};
[rxGroundHeightWrtTxInM, rxHeightWrtTxInM] = deal(cell(numTracks,1));
for idxTrack = 1:numTracks
    % RX location.
    fullPathRxLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['rxLoc_track_', num2str(idxTrack), '.csv']);
    curHeaderCell = {'utmX', 'utmY', 'utmZone', ...
        'lat', 'lon', 'altFromUsgs', 'rxHeightWrtTXInM'};
    curXsRx = pathLossUtmXYHs{idxTrack}(:,1);
    curYsRx = pathLossUtmXYHs{idxTrack}(:,2);
    curAltsUsgsRx = getUsgsAltInM(curXsRx, curYsRx);
    curRxHeightWrtTxInM = (RX_HEIGHT_M+curAltsUsgsRx)...
        -(TX_ALT_USGS+TX_HEIGHT_M);
    
    rxGroundHeightWrtTxInM{idxTrack} ...
        = curAltsUsgsRx-(TX_ALT_USGS+TX_HEIGHT_M);
    rxHeightWrtTxInM{idxTrack} = curRxHeightWrtTxInM;
    
    curData = [num2cell([curXsRx, ...
        curYsRx]), ...
        cellstr(pathLossUtmZones{idxTrack}), ...
        ...
        num2cell([contiPathLossesWithGpsInfo{idxTrack}(:,2), ...
        contiPathLossesWithGpsInfo{idxTrack}(:,3), ...
        curAltsUsgsRx, ...
        curRxHeightWrtTxInM])];
    writeToCsvWithHeader(fullPathRxLocs, curHeaderCell, ...
        curData, dataStrFormatter);
    
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
    writeToCsvWithHeader(fullPathRxMeas, curHeaderCell, ...
        curData, dataStrFormatter);
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

curFigTitle = 'Simulation Overview in UTM';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'simOverviewUtmXY.png');

curFig = figure; hold on;
hGround = plot3(VEG_AREA_IMG_META.XS(:), VEG_AREA_IMG_META.YS(:), ...
    groundHeightWrtTXInM(:), '.', 'MarkerSize', 1, 'Color', lightGray);
hTx = plot3(xTx,yTx,0, 'gx');
for idxTrack = 1:numTracks
    hRx = plot3(pathLossUtmXYHs{idxTrack}(:,1), ...
        pathLossUtmXYHs{idxTrack}(:,2), ...
        rxHeightWrtTxInM{idxTrack}, ...
        'r.', 'MarkerSize', 5);
end
hTrees = plot3(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), ...
    trunkGroundHeightWrtTXInM+trunkTreeHeightInM, '.b');
grid on; view(2); axis tight; axis equal;
xlabel('utmX'); ylabel('utmY'); uistack(hGround, 'bottom');
legend([hTx, hRx, hTrees, hGround], ...
    'TX', 'RX', 'Trees', 'Ground', 'Location', 'SouthEast');

title(curFigTitle)
saveas(curFig, curFigPath);

curFigTitle = 'Simulation Overview in 3D';
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'simOverviewUtmXYH.png');

curFig = figure; hold on;
downSampleRate = 30;
fctDownSample = @(x) x(1:downSampleRate:end, 1:downSampleRate:end);
vegXsLess = fctDownSample(VEG_AREA_IMG_META.XS);
vegYsLess = fctDownSample(VEG_AREA_IMG_META.YS);
vegGroundHsWrtTxLess = fctDownSample(groundHeightWrtTXInM);
hGround = plot3(vegXsLess(:), vegYsLess(:), ...
    vegGroundHsWrtTxLess(:), '.', 'MarkerSize', 1, 'Color', lightGray);
plot3([xTx xTx],[yTx yTx],[-TX_HEIGHT_M 0], 'g-');
hTx = plot3(xTx,yTx,0, 'gx');
for idxTrack = 1:numTracks
    plot3([pathLossUtmXYHs{idxTrack}(:,1) ...
        pathLossUtmXYHs{idxTrack}(:,1)]', ...
        [pathLossUtmXYHs{idxTrack}(:,2) ...
        pathLossUtmXYHs{idxTrack}(:,2)]', ...
        [rxHeightWrtTxInM{idxTrack} rxGroundHeightWrtTxInM{idxTrack}]', ...
        'r-');
    hRx = plot3(pathLossUtmXYHs{idxTrack}(:,1), ...
        pathLossUtmXYHs{idxTrack}(:,2), ...
        rxHeightWrtTxInM{idxTrack}, ...
        'r.', 'MarkerSize', 5);
end
plot3([treeUtmXYsExt(:,1) treeUtmXYsExt(:,1)]', ...
    [treeUtmXYsExt(:,2) treeUtmXYsExt(:,2)]', ...
    [trunkGroundHeightWrtTXInM ...
    trunkGroundHeightWrtTXInM+trunkTreeHeightInM]', '-b');
hTrees = plot3(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), ...
    trunkGroundHeightWrtTXInM+trunkTreeHeightInM, 'ob', 'MarkerSize', 3);
grid on; axis tight; axis equal; view(3);
xlabel('utmX'); ylabel('utmY'); ylabel('Height');
legend([hTx, hRx, hTrees, hGround], ...
    'TX', 'RX', 'Trees', 'Ground', 'Location', 'SouthEast');

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
    hTrunks = plot3(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2), ...
        ones(length(treeUtmXYsExt(:,1)),1) ...
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
    hTrunks = plot(treeUtmXYsExt(:,1), treeUtmXYsExt(:,2),...
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
    plot(markLocs(indicesForTrunksWithNegHeight,2), ...
        markLocs(indicesForTrunksWithNegHeight,1), 'rx');
    plot_google_map('MapType', 'satellite');
    xticks([]); yticks([]);
    title('Trunk Locations with Negtive Tree Height');
end

%% Figures According to the Exported .csv Files

% foliageAreaGrid.csv

% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "utmX", "utmY", "Var5", "Var6", ...
    "Var7", "Var8", "groundHeightWrtTXInM", "Var10", "Var11"];
opts.SelectedVariableNames = ["utmX", "utmY", "groundHeightWrtTXInM"];
opts.VariableTypes = ["string", "string", "double", "double", "string", ...
    "string", "string", "string", "double", "string", "string"];
opts = setvaropts(opts, [1, 2, 5, 6, 7, 8, 10, 11], ...
    "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 5, 6, 7, 8, 10, 11], ...
    "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
foliageAreaGridExp = readtable(fullPathFoliageAreaCsv, opts);
vegXsExp = foliageAreaGridExp.utmX;
vegYsExp = foliageAreaGridExp.utmY;
vegGroundHsWrtTxExp = foliageAreaGridExp.groundHeightWrtTXInM;

% Clear temporary variables
clear opts

% txLoc.csv
txLocExp = readtable(fullPathTxLoc);
xTxExp = txLocExp.utmX;
yTxExp = txLocExp.utmY;
txHExp = txLocExp.heightAboveGroundInM;

% trunkLocs.csv
trunkLocsExp = readtable(fullPathTrunkLocCsv);
treeUtmXYsExtExp = [trunkLocsExp.utmX trunkLocsExp.utmY];
trunkGroundHeightWrtTXExp = trunkLocsExp.groundHeightWrtTXInM;
trunkTreeHeightExp = trunkLocsExp.treeHeightInM;

% Overview figure.
curFigTitle = {'Simulation Overview in 3D'; 'Based on Exported Files'};
curFigPath = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'simOverviewUtmXYHByCsv.png');

downSampleRate = 30;
fctDownSample = @(x) x(1:downSampleRate:end, 1:downSampleRate:end);
vegXsLessExp = fctDownSample(vegXsExp);
vegYsLessExp = fctDownSample(vegYsExp);
vegGroundHsWrtTxLessExp = fctDownSample(vegGroundHsWrtTxExp);

curFig = figure; hold on;
hGround = plot3(vegXsLessExp(:), vegYsLessExp(:), ...
    vegGroundHsWrtTxLessExp(:), '.', 'MarkerSize', 1, 'Color', lightGray);
plot3([xTxExp xTxExp],[yTxExp yTxExp],[-txHExp 0], 'g-');
hTx = plot3(xTxExp,yTxExp,0, 'gx');
for idxTrack = 1:numTracks
    % rxLoc_track_#.csv
    fullPathRxLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['rxLoc_track_', num2str(idxTrack), '.csv']);
    rxLocExp = readtable(fullPathRxLocs);
    pathLossUtmXYsExp = [rxLocExp.utmX, rxLocExp.utmY];
    rxHeightWrtTxExp = rxLocExp.rxHeightWrtTXInM;
    
    plot3([pathLossUtmXYsExp(:,1) ...
        pathLossUtmXYsExp(:,1)]', ...
        [pathLossUtmXYsExp(:,2) ...
        pathLossUtmXYsExp(:,2)]', ...
        [rxHeightWrtTxExp rxGroundHeightWrtTxInM{idxTrack}]', ...
        'r-');
    hRx = plot3(pathLossUtmXYsExp(:,1), ...
        pathLossUtmXYsExp(:,2), ...
        rxHeightWrtTxExp, ...
        'r.', 'MarkerSize', 5);
end
plot3([treeUtmXYsExtExp(:,1) treeUtmXYsExtExp(:,1)]', ...
    [treeUtmXYsExtExp(:,2) treeUtmXYsExtExp(:,2)]', ...
    [trunkGroundHeightWrtTXExp ...
    trunkGroundHeightWrtTXExp+trunkTreeHeightExp]', '-b');
hTrees = plot3(treeUtmXYsExtExp(:,1), treeUtmXYsExtExp(:,2), ...
    trunkGroundHeightWrtTXExp+trunkTreeHeightExp, 'ob', 'MarkerSize', 3);
grid on; axis tight; axis equal; view(3);
xlabel('utmX'); ylabel('utmY'); ylabel('Height');
legend([hTx, hRx, hTrees, hGround], ...
    'TX', 'RX', 'Trees', 'Ground', 'Location', 'SouthEast');

title(curFigTitle)
saveas(curFig, curFigPath);

%% Close Figures
close all;

% EOF