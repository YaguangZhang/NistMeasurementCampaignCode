% ESTIMATEFOLIAGEATTENUATION Estimate the extra path losses caused by trees
% for the NIST dataset.
%
% Yaguang Zhang, Purdue, 04/16/2018

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
    'PostProcessingResults', 'FoliageAttenuationEstimation');
ABS_PATH_TO_TREE_LOCS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'Data', 'NIST foliage analysis tree locations.csv');
% The loaded tree info only has lat and lon. We will fetch alt accordingly
% and save the results here.
pathToSaveTreeLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'treeLocs.mat');
pathToSaveUtmInfo = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'utmInfoForPathLossesAndTrees.mat');
pathToSaveNumsOfTreesInFirstFresnel = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'numsTreesInFirstFresnel.mat');
pathToSaveFreeSpacePathLosses = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'freeSpacePathLossResults.mat');

% Reuse results from evalPathLossesForContiTracks.m and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

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

assert(exist(ABS_PATH_TO_PATH_LOSSES, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/3_PathLossComputation/evalPathLossesForContiTracks.m first.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');

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

disp('    Done!')

%% Load the Tree Info

disp(' ')
disp('    Loading tree location information...')

treeLocations = loadGpsMarkersWithAlt(ABS_PATH_TO_TREE_LOCS, 'Marker*', ...
    pathToSaveTreeLocs);
% Keep trying to fetch alts from Google until all needed info is available.
numRetrial = 0;
while any(isnan(treeLocations(:,3)))
    numRetrial = numRetrial+1;
    disp(['Not all alts retrieved. Retry again ... (#', ...
        num2str(numRetrial), ')']);
    treeLocations = loadGpsMarkersWithAlt(ABS_PATH_TO_TREE_LOCS, 'Marker*', ...
        pathToSaveTreeLocs);
end
disp('    Done!')

%% Convert All GPS Samples to UTM

disp(' ')
disp('    Converting GPS samples to UTM ...')

% Path loss samples contiPathLossesWithGpsInfo.
numOfSeries = length(contiPathLossesWithGpsInfo);
[pathLossUtmXYHs,pathLossUtmZones] ...
    = deal(cell(numOfSeries,1));
for idxS = 1:numOfSeries
    [xs,ys,zones] ...
        = deg2utm(contiPathLossesWithGpsInfo{idxS}(:,2), ...
        contiPathLossesWithGpsInfo{idxS}(:,3));
    pathLossUtmXYHs{idxS} = [xs ys contiPathLossesWithGpsInfo{idxS}(:,4)];
    pathLossUtmZones{idxS} = zones;
end
% Tree locations.
[xs,ys,treeUtmZones] = deg2utm(treeLocations(:,1), ...
    treeLocations(:,2));
treeUtmXYHs = [xs ys treeLocations(:,3)];

% Tx location.
[xTx, yTx, txUtmZone] = deg2utm(TX_LAT, TX_LON);

% We are not yet able to deal with points from different UTM zones.
allRelevantZones ...
    = num2cell([vertcat(pathLossUtmZones{:});treeUtmZones], 2);
assert(all(strcmp(allRelevantZones, txUtmZone)), ...
    'We are not yet able to deal with points from different UTM zones!');

save(pathToSaveUtmInfo, ...
    'xTx', 'yTx', 'txUtmZone', 'treeLocations', ...
    'pathLossUtmXYHs', 'pathLossUtmZones', ...
    'treeUtmXYHs', 'treeUtmZones');

disp('    Done!')

%% Plot an Overview in UTM

disp(' ')
disp('    Generating overview in UTM ...')

pathToSaveOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'OverviewInUtm');
hOverviewInUtm = figure;
hold on;
allPathLossUtmXYHs = vertcat(pathLossUtmXYHs{:});
allPathLosses = vertcat(contiPathLossesWithGpsInfo{:});
allPathLosses = allPathLosses(:,1);
boolsFinitePathLosses = ~isinf(allPathLosses);
% Trees.
hTrees = plot3(treeUtmXYHs(:,1), treeUtmXYHs(:,2), treeUtmXYHs(:,3), 'g*');
% Pathlosses.
plot3k(allPathLossUtmXYHs(boolsFinitePathLosses,:), ...
    'ColorData', allPathLosses(boolsFinitePathLosses));
title('Samples in 3D (Colored by path losses) with trees');
xlabel('x (m)'); ylabel('y (m)'); zlabel('altitude (m)');

saveas(hOverviewInUtm, [pathToSaveOverviewPlot, '.fig']);
saveas(hOverviewInUtm, [pathToSaveOverviewPlot, '.png']);

disp('    Done!')

%% Count the Number of Trees in the Fresii for Each Path Loss Measurement
% We will also append the results into .mat file for future use.

disp(' ')
disp('    Counting trees in Fresnel zones ...')

numsOfTreesInFirstFresnel = cell(numOfSeries, 1);
tx3D = [xTx, yTx, TX_ALT];
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);
    [curNumSamps, ~] = size(pathLossUtmXYHs{idxS});
    numsOfTreesInFirstFresnel{idxS} = nan(curNumSamps, 1);
    
    for idxSamp = 1:curNumSamps
        rx3D = pathLossUtmXYHs{idxS}(idxSamp, :);
        numsOfTreesInFirstFresnel{idxS}(idxSamp) ...
            = countNumOfTreesInFirstFresnelZone(tx3D, rx3D, ...
            treeUtmXYHs, F_C_IN_GHZ);
    end
end

save(pathToSaveNumsOfTreesInFirstFresnel, 'numsOfTreesInFirstFresnel');

disp('    Done!')

%% Compute Free-Space Path Losses

disp(' ')
disp('    Computing free-space path losses ...')

% TX-RX distances.
txRxLosDists = cell(numOfSeries, 1);
% Free-space path losses in dB.
freeSpacePathLosses = cell(numOfSeries, 1);
% The excessive path losses in dB (measured path losses -  free space path
% losses).
exceLossRefFreeSpace = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);
    [curNumSamps, ~] = size(pathLossUtmXYHs{idxS});
    [freeSpacePathLosses{idxS}, txRxLosDists{idxS}] ...
        = computeFreeSpacePathLosses(tx3D, ...
        pathLossUtmXYHs{idxS}, F_C_IN_GHZ);
    exceLossRefFreeSpace{idxS} = contiPathLossesWithGpsInfo{idxS}(:,1) ...
        - freeSpacePathLosses{idxS};
end

save(pathToSaveFreeSpacePathLosses, 'freeSpacePathLosses', ...
    'exceLossRefFreeSpace');

disp('    Done!')

%% Plot Figures for the Tree Numbers

disp(' ')
disp('    Generating plots for tree numbers ...')

pathToSaveTreeNumPlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'TreeNums');

% Tree numbers on a 2D map.
hTreeNumOnMap = figure;
hold on;
allTreeNumsInFirstFresnel = vertcat(numsOfTreesInFirstFresnel{:});
allContiPathLossesWithGpsInfo = vertcat(contiPathLossesWithGpsInfo{:});
% We need (lon, lat) for plotting.
allRxGpsLocs = allContiPathLossesWithGpsInfo(:, [3 2]);
% TX.
hTx = plot3(TX_LON, TX_LAT, 0, 'kx');
% Trees.
hTrees = plot3(treeLocations(:,2), treeLocations(:,1), ...
    treeLocations(:,3).*0, 'g*');
% Tree numbers in 1st Fresnel zone for each RX location.
plot3k([allRxGpsLocs allTreeNumsInFirstFresnel], 'PlotType', 'stem');
plot_google_map('MapType', 'satellite');
view(95, 70);
title({'Number of Trees in the 1st Fresnel Zone', 'for Each RX Location'});
xlabel(' '); ylabel(' '); xticks([]); yticks([]);
zlabel('Number of Trees'); 

saveas(hTreeNumOnMap, [pathToSaveTreeNumPlots, 'OnMap.fig']);
saveas(hTreeNumOnMap, [pathToSaveTreeNumPlots, 'OnMap.png']);

% Tree numbers and excessive path losses over TX-RX distance.
hTreeNumAndExceLoss = figure;
hold on;
allTxRxLosDists = vertcat(txRxLosDists{:});
[sortedLosDist, indicesForSortedDist] = sort(allTxRxLosDists); 
allExceLossRefFreeSpace = vertcat(exceLossRefFreeSpace{:});
% Tree numbers.
yyaxis left;
plot(sortedLosDist, allTreeNumsInFirstFresnel(indicesForSortedDist));
ylabel('Number of Trees');
% Excessive path losses.
yyaxis right;
plot(sortedLosDist, allExceLossRefFreeSpace(indicesForSortedDist));
title({'Number of Trees in the 1st Fresnel Zone', ...
    'with Excessive Path Losses (Over Free-Space Losses)'});
xlabel('RX Distance (m)'); ylabel('Excessive Path Loss (dB)'); grid on;

saveas(hTreeNumAndExceLoss, ...
    [pathToSaveTreeNumPlots, 'AndExcePathLosses.fig']);
saveas(hTreeNumAndExceLoss, ...
    [pathToSaveTreeNumPlots, 'AndExcePathLosses.png']);

% Measured path losses and free-space path losses.
hMeasAndFreeSpaceLosses = figure;
allFreeSpacePathLosses = vertcat(freeSpacePathLosses{:});
% Measured path losses.
yyaxis left;
semilogx(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1));
ylabel('Measured Path Loss (dB)');
% Free-space path losses
yyaxis right;
semilogx(sortedLosDist, allFreeSpacePathLosses(indicesForSortedDist));
title({'Measured and Free-Space Path Losses'});
xlabel('RX Distance (m)'); ylabel('Free-Space Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLosses, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.fig'));
saveas(hMeasAndFreeSpaceLosses, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.png'));

% Measured path losses and free-space path losses sharing the same y axis.
hMeasAndFreeSpaceLossesSameY = figure;
hold on;
% Measured path losses.
hMeas = plot(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1));
% Free-space path losses
hExce = plot(sortedLosDist, allFreeSpacePathLosses(indicesForSortedDist));
set(gca, 'XScale', 'log')
legend([hMeas, hExce], 'Measured', 'Free-Space');
title({'Measured and Free-Space Path Losses'});
xlabel('RX Distance (m)'); ylabel('Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLossesSameY, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.fig'));
saveas(hMeasAndFreeSpaceLossesSameY, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.png'));

% Excessive path loss over number of trees in the 1st Fresnel zone.
hExcePathLossOverTreeNum = figure;
plot(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, '.');
title({'Excessive Path Losses vs. Number of Trees in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excessive Path Losses (dB)'); 
axis tight; grid on;

saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.fig'));
saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.png'));

% Excessive path loss over log number of trees in the 1st Fresnel zone.
hExcePathLossOverLogTreeNum = figure;
semilogx(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, '.');
title({'Excessive Path Losses vs. Log Number of Trees in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excessive Path Losses (dB)'); 
axis tight; grid on;

saveas(hExcePathLossOverLogTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverLogTreeNum.fig'));
saveas(hExcePathLossOverLogTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverLogTreeNum.png'));

disp('    Done!')

% EOF