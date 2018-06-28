% ESTIMATEFOLIAGEATTENUATION Estimate the extra path losses caused by trees
% for the NIST dataset.
%
%   This script is based on:
%
%       4_FoliageAttenuationEstimation/estimateFoliageAttenuation.m
%
%   Tree locations manually marked on a map with satellite images and lidar
%   data for NIST, instead of the records made by the GPS marker Android
%   app is used here.
%
% Yaguang Zhang, Purdue, 05/31/2018

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Need functions from 4_FoliageAttenuationEstimation.
addpath(fullfile(pwd, '4_FoliageAttenuationEstimation'));

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'FoliageAttenuationEstimation_ManualTreeLocs');
ABS_PATH_TO_MANUAL_TREE_LOCS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'ManuallyLocateTrees', 'treeLocs.mat');
% The loaded tree info only has lat and lon. We will fetch alt accordingly
% and save the results here.
pathToSaveTreeLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'treeLocsWithAlts.mat');
pathToSaveUtmInfo = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'utmInfoForPathLossesAndTrees.mat');
pathToSaveFoliageAttenResults = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'foliageAttenAnalysisResults.mat');

% Reuse results from evalPathLossesForContiTracks.m and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

GOOGLE_MAPS_API = 'AIzaSyDlkaE_QJxvRJpTutWG0N-LCvoT0e7FPHE';

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
disp('      - treeLocs.mat')
disp('      - contiPathLossesWithGpsInfo.mat')
disp('      - txInfoLogs.mat')

assert(exist(ABS_PATH_TO_MANUAL_TREE_LOCS, 'file')==2, ...
    'Couldn''t find treeLocs.mat! Please run PostProcessing/7_ManuallyLocateTrees/manuallyLocateTrees.m first.');
assert(exist(ABS_PATH_TO_PATH_LOSSES, 'file')==2, ...
    'Couldn''t find contiPathLossesWithGpsInfo.mat! Please run PostProcessing/3_PathLossComputation/evalPathLossesForContiTracks.m first.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'markLocs'.
load(ABS_PATH_TO_MANUAL_TREE_LOCS);
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
disp('    Fetching altitudes for tree locations ...')

[numTrees, ~] = size(markLocs);
treeLocations = [markLocs nan(numTrees,1)];
% Keep trying to fetch alts from Google until all needed info is available.
numRetrial = 0;
curBoolsNanAlt = isnan(treeLocations(:,3));
while any(curBoolsNanAlt)
    numRetrial = numRetrial+1;
    disp(['Not all alts retrieved. Retry again ... (#', ...
        num2str(numRetrial), ')']);

    treeLocations(curBoolsNanAlt, 3) = getElevationsFromGoogle( ...
        treeLocations(curBoolsNanAlt, 1), ...
        treeLocations(curBoolsNanAlt, 2), GOOGLE_MAPS_API);
    curBoolsNanAlt = isnan(treeLocations(:,3));
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

%% Plot an Overview in UTM for the Path Losses

disp(' ')
disp('    Generating overview in UTM for path losses ...')

pathToSaveOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'OverviewInUtm');
hOverviewInUtm = figure;
hold on;
allPathLossUtmXYHs = vertcat(pathLossUtmXYHs{:});
allPathLosses = vertcat(contiPathLossesWithGpsInfo{:});
allPathLosses = allPathLosses(:,1);
boolsFinitePathLosses = ~isinf(allPathLosses);
% Trees.
hTrees = plot3(treeUtmXYHs(:,1), treeUtmXYHs(:,2), ...
    treeUtmXYHs(:,3).*min(allPathLosses), 'g*');
% Pathlosses.
plot3k(allPathLossUtmXYHs(boolsFinitePathLosses,:), ...
    'ColorData', allPathLosses(boolsFinitePathLosses));
title('Samples in 3D (Colored by Path Losses) with Trees');
view(0, 90);
xlabel('x (m)'); ylabel('y (m)'); zlabel('Path loss (dB)');

saveas(hOverviewInUtm, [pathToSaveOverviewPlot, '.fig']);
saveas(hOverviewInUtm, [pathToSaveOverviewPlot, '.png']);

disp('    Done!')

%% Plot an Overview in UTM for GPS Samples

disp(' ')
disp('    Generating overview in UTM for GPS samples ...')

pathToSaveGpsOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'OverviewInUtmGps');
hOverviewInUtmGps = figure;
hold on;
% Trees.
hTrees = plot3(treeUtmXYHs(:,1), treeUtmXYHs(:,2), treeUtmXYHs(:,3), 'g*');
% TX.
hTx = plot3(xTx, yTx, TX_ALT+TX_HEIGHT_M, 'kx');
% RX locations.
allPathLossUtmXYHsShifted = allPathLossUtmXYHs;
allPathLossUtmXYHsShifted(:,3) ...
    = allPathLossUtmXYHsShifted(:,3)+RX_HEIGHT_M;
hRxs = plot3k(allPathLossUtmXYHsShifted);
title('Samples in 3D (Colored by Altitudes) with Trees');
xlabel('x (m)'); ylabel('y (m)'); zlabel('altitude (m)');

saveas(hOverviewInUtmGps, [pathToSaveGpsOverviewPlot, '.fig']);
saveas(hOverviewInUtmGps, [pathToSaveGpsOverviewPlot, '.png']);

disp('    Done!')

%% Count the Number of Trees in the Fresii for Each Path Loss Measurement
% We will also append the results into .mat file for future use.

disp(' ')
disp('    Counting trees in Fresnel zones ...')

numsOfTreesInFirstFresnel = cell(numOfSeries, 1);
tx3D = [xTx, yTx, TX_ALT+TX_HEIGHT_M];
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);
    [curNumSamps, ~] = size(pathLossUtmXYHs{idxS});
    numsOfTreesInFirstFresnel{idxS} = nan(curNumSamps, 1);

    for idxSamp = 1:curNumSamps
        rx3D = pathLossUtmXYHs{idxS}(idxSamp, :);
        % Considering the RX height over the ground.
        rx3D(:, 3) = rx3D(:, 3)+RX_HEIGHT_M;
        numsOfTreesInFirstFresnel{idxS}(idxSamp) ...
            = countNumOfTreesInFirstFresnelZone(tx3D, rx3D, ...
            treeUtmXYHs, F_C_IN_GHZ);
    end
end

save(pathToSaveFoliageAttenResults, 'numsOfTreesInFirstFresnel');

disp('    Done!')

%% Compute Free-Space Path Losses

disp(' ')
disp('    Computing free-space path losses ...')

% TX-RX distances.
txRxLosDists = cell(numOfSeries, 1);
% Free-space path losses in dB.
freeSpacePathLosses = cell(numOfSeries, 1);
% The excess path losses in dB (measured path losses -  free space path
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

save(pathToSaveFoliageAttenResults, 'freeSpacePathLosses', ...
    'exceLossRefFreeSpace', '-append');

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
% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hTreeNumOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

saveas(hTreeNumOnMap, [pathToSaveTreeNumPlots, 'OnMap.fig']);
saveas(hTreeNumOnMap, [pathToSaveTreeNumPlots, 'OnMap.png']);

% Tree numbers and excess path losses over TX-RX distance.
hTreeNumAndExceLoss = figure;
hold on;
allTxRxLosDists = vertcat(txRxLosDists{:});
[sortedLosDist, indicesForSortedDist] = sort(allTxRxLosDists);
allExceLossRefFreeSpace = vertcat(exceLossRefFreeSpace{:});
% Tree numbers.
yyaxis left;
plot(sortedLosDist, allTreeNumsInFirstFresnel(indicesForSortedDist));
ylabel('Number of Trees');
% Excess path losses.
yyaxis right;
plot(sortedLosDist, allExceLossRefFreeSpace(indicesForSortedDist));
title({'Number of Trees in the 1st Fresnel Zone', ...
    'with Excess Path Losses (Over Free-Space Losses)'});
xlabel('RX Distance (m)'); ylabel('Excess Path Loss (dB)'); grid on;

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

% Excess path loss over number of trees in the 1st Fresnel zone.
hExcePathLossOverTreeNum = figure;
plot(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, '.');
title({'Excess Path Losses vs. Number of Trees in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excess Path Losses (dB)');
axis tight; grid on;

saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.fig'));
saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.png'));

% Excess path loss over log number of trees in the 1st Fresnel zone.
hExcePathLossOverLogTreeNum = figure;
semilogx(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, '.');
title({'Excess Path Losses vs. Log Number of Trees in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excess Path Losses (dB)');
axis tight; grid on;

saveas(hExcePathLossOverLogTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverLogTreeNum.fig'));
saveas(hExcePathLossOverLogTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverLogTreeNum.png'));

disp('    Done!')

%% Plot an Overview Comparing GPS Alts with Google Elevations
disp(' ')
disp('    Generating overview in UTM for comparing GPS alts with Google eles ...')

if contiPathLossesExtraInfo.FLAG_USE_GOOGLE_FOR_ALT

    allContiOriAltsShifted ...
        = vertcat(contiPathLossesExtraInfo.contiOriAlts{:}) + RX_HEIGHT_M;

    pathToSaveGpsOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        'OverviewInUtmGpsComp');
    hOverviewInUtmGpsComp = figure;
    hold on;
    % Trees.
    hTrees = plot3(treeUtmXYHs(:,1), treeUtmXYHs(:,2), treeUtmXYHs(:,3), 'g*');
    % TX.
    hTx = plot3(xTx, yTx, TX_ALT+TX_HEIGHT_M, 'kx');
    % RX locations according to Google.
    hRxsG = plot3(allPathLossUtmXYHsShifted(:,1), ...
        allPathLossUtmXYHsShifted(:,2), allPathLossUtmXYHsShifted(:,3), ...
        'b.');
    % RX locations according to GPS samples.
    hRxs = plot3(...
        allPathLossUtmXYHsShifted(:,1), allPathLossUtmXYHsShifted(:,2), ...
        allContiOriAltsShifted, 'r.');
    title('Comparing GPS Alts with Google Elevations in 3D with Trees');
    legend([hTx, hTrees, hRxs, hRxsG], 'TX', 'Trees', ...
        'RX with GPS Alts', 'RX with Google Eles');
    transparentizeCurLegends; grid minor; view(40, 5);
    xlabel('UTM x (m)'); ylabel('UTM y (m)'); zlabel('Altitude (m)');

    saveas(hOverviewInUtmGpsComp, [pathToSaveGpsOverviewPlot, '.fig']);
    saveas(hOverviewInUtmGpsComp, [pathToSaveGpsOverviewPlot, '.png']);

    % Plot the difference (alts - eles) vs distance.
    pathToSaveGpsAltAndGoogleEleDiffVsDist = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        'OverviewGpsAltAndGoogleEleDiffVsDist');
    hOverviewGpsAltAndGoogleEleDiffVsDist = figure;
    plot(sortedLosDist, allContiOriAltsShifted(indicesForSortedDist) ...
        - allPathLossUtmXYHsShifted(indicesForSortedDist,3));
    title('Difference between GPS Alts and Google Elevations over Distance');
    xlabel('Distance (m)'); ylabel('GPS Alts - Google Eles (m)');
    grid minor; axis tight;

    saveas(hOverviewGpsAltAndGoogleEleDiffVsDist, ...
        [pathToSaveGpsAltAndGoogleEleDiffVsDist, '.fig']);
    saveas(hOverviewGpsAltAndGoogleEleDiffVsDist, ...
        [pathToSaveGpsAltAndGoogleEleDiffVsDist, '.png']);

    % Plot the difference (alts - eles).
    pathToSaveGpsAltAndGoogleEleDiffBin = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        'OverviewGpsAltAndGoogleEleDiffBin');
    hOverviewGpsAltAndGoogleEleDiffBin = figure;
    histogram( allContiOriAltsShifted - allPathLossUtmXYHsShifted(:,3));

    title('Histogram for Difference between GPS Alts and Google Elevations');
    xlabel('GPS Alts - Google Eles (m)'); ylabel('#'); grid minor;

    saveas(hOverviewGpsAltAndGoogleEleDiffBin, ...
        [pathToSaveGpsAltAndGoogleEleDiffBin, '.fig']);
    saveas(hOverviewGpsAltAndGoogleEleDiffBin, ...
        [pathToSaveGpsAltAndGoogleEleDiffBin, '.png']);

else
    disp('        Google elevations were not retreived. Abortting...')
end

disp('    Done!')

%% Plot Path Loss Results with Antenna Gains.
disp(' ')
disp('    Plotting path loss results with antenna gains ...')

alphaForResultsExcludingAntGains = 0.8;

allContiPathTxGains = vertcat(contiPathLossesExtraInfo.contiTxGains{:});
allContiPathRxGains = vertcat(contiPathLossesExtraInfo.contiRxGains{:});

pathToSavePathLossWithAntGains = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'PathLossWithAntGainsOverDist');
hPathLossWithAntGainsOverDist = figure; hold on;
hMeas = plot(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1), '--');
hMeasMinusRxGains = plot(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1) ...
    - allContiPathRxGains(indicesForSortedDist), '-');
hMeasMinusRxAndTxGains = plot(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1) ...
    - allContiPathRxGains(indicesForSortedDist) ...
    - allContiPathTxGains(indicesForSortedDist), '-.');
hMeasMinusRxGains.Color(4) = alphaForResultsExcludingAntGains;
hMeasMinusRxAndTxGains.Color(4) = alphaForResultsExcludingAntGains;
legend([hMeas, hMeasMinusRxGains, hMeasMinusRxAndTxGains], ...
    'Measurements', 'Measurements excluding RX gain', ...
    'Measurements excluding RX and TX gain');
transparentizeCurLegends;
title({'Path Losses with Antenna Gains over Distance'});
xlabel('Distance (m)'); ylabel('Path Losses (dB)');
axis tight; grid minor;

saveas(hPathLossWithAntGainsOverDist, ...
    [pathToSavePathLossWithAntGains, '.fig']);
saveas(hPathLossWithAntGainsOverDist, ...
    [pathToSavePathLossWithAntGains, '.png']);

% The excess path loss vs tree number for each track.
numTracks = length(exceLossRefFreeSpace);
for idxTrack = 1:numTracks
    curExceLossRefFreeSpace = exceLossRefFreeSpace{idxTrack};
    curNumsOfTreesInFirstFresnel = numsOfTreesInFirstFresnel{idxTrack};

    hExcePathLossOverTreeNumCurTrack = figure; hold on;
    curNumPtsToShow = length(curNumsOfTreesInFirstFresnel);
    scatter(curNumsOfTreesInFirstFresnel, curExceLossRefFreeSpace, ...
        ones(curNumPtsToShow, 1) .* 8, ...
        'MarkerFaceColor','b','MarkerEdgeColor','none',...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    title({'Excess Path Losses vs. Number of Trees in the 1st Fresnel Zone', ...
        ['Track ', num2str(idxTrack)]});
    xlabel('Number of Trees'); ylabel('Excess Path Losses (dB)');
    axis tight; grid minor;
    % Extend the x range a little bit to better show the dots on x = 0.
    curAxis = axis; axis([-0.5 curAxis(2:end)]);

    saveas(hExcePathLossOverTreeNumCurTrack, ...
        fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['ExceLossOverTreeNumTrack_', num2str(idxTrack), '.fig']));
    saveas(hExcePathLossOverTreeNumCurTrack, ...
        fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['ExceLossOverTreeNumTrack_', num2str(idxTrack), '.png']));
end

disp('    Done!')

%% Apply the Partition-Dependent Model
% Essentially we will get one fixed value for the extra path loss caused by
% each tree inside the 1st Fresnel zone.

disp(' ')
disp('    Computing the extra path loss caused by each tree in the 1st Fresnel zone...')

excePathLossPerTree ...
    = pinv(allTreeNumsInFirstFresnel)*allExceLossRefFreeSpace;

disp('    Making predictions accordingly ...')
shiftedFreeSpacePathLosses = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);

    shiftedFreeSpacePathLosses{idxS} = freeSpacePathLosses{idxS} ...
        + numsOfTreesInFirstFresnel{idxS}.*excePathLossPerTree;
end

save(pathToSaveFoliageAttenResults, 'excePathLossPerTree', ...
    'shiftedFreeSpacePathLosses', '-append');

%% Plot the Pathloss Results

allShiftedFreeSpacePathLosses = vertcat(shiftedFreeSpacePathLosses{:});

% Measured path losses and free-space path losses over linear distance.
hMeasAndFreeSpaceLossesLinearDist = figure; hold on;
allFreeSpacePathLosses = vertcat(freeSpacePathLosses{:});
% Measured path losses.
hMeas = plot(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1), '-');
% Free-space path losses.
hFreeSpace = plot(sortedLosDist, ...
    allFreeSpacePathLosses(indicesForSortedDist), '-');
% Shifted free-space path losses.
hShiftedFreeSpace ...
    = plot(sortedLosDist, ...
    allShiftedFreeSpacePathLosses(indicesForSortedDist), '-');
legend([hMeas, hFreeSpace, hShiftedFreeSpace], ...
    'Measurements', 'Free Space Path Loss', ...
    'Shifted Free Space Path Loss', 'Location', 'southeast');
axis tight; transparentizeCurLegends;

% RMSE
rmseShiftVsMeas = sqrt(sum((allShiftedFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allShiftedFreeSpacePathLosses));

title(['Foilage Analysis Results: RMSE = ', ...
    num2str(rmseShiftVsMeas, '%.4f'), ' dB']);
xlabel('RX Distance (m)'); ylabel('Free-Space Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLossesLinearDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLinearDist.fig'));
saveas(hMeasAndFreeSpaceLossesLinearDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLinearDist.png'));

% Excess path loss over number of trees in the 1st Fresnel zone with
% some analysis information overlaid.
hExcePathLossOverTreeNumMore = figure; hold on;
numPtsToShow = length(allTreeNumsInFirstFresnel);
hExcePathLosses ...
    = scatter(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, ...
    ones(numPtsToShow, 1) .* 8, ...
    'MarkerFaceColor','b','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% Also show the median, mean and fitted lines.
treeNumValues = min(allTreeNumsInFirstFresnel) ...
    :max(allTreeNumsInFirstFresnel);
excePathLossMean = arrayfun(@(n) ...
    mean(allExceLossRefFreeSpace(allTreeNumsInFirstFresnel==n)), ...
    treeNumValues);
excePathLossMedian = arrayfun(@(n) ...
    median(allExceLossRefFreeSpace(allTreeNumsInFirstFresnel==n)), ...
    treeNumValues);
fittedCurve = fit(allTreeNumsInFirstFresnel, allExceLossRefFreeSpace, ...
    'smoothingspline');
% The fitted line will be stored as [x(1), y(1); x(2) y(2)].
fittedLine = [[0; treeNumValues(end)], ...
    [0; treeNumValues(end)].*excePathLossPerTree];
% Also plot a flat line (a constant offset with slope 0).
meanOfAllExcePathLosses = mean(allExceLossRefFreeSpace);
% Plots.
curLineWidth = 1;
hExcePathLossMeans = plot(treeNumValues, excePathLossMean, 'xr-.', ...
    'LineWidth', curLineWidth);
hExcePathLossMedians = plot(treeNumValues, excePathLossMedian, '+b--', ...
    'LineWidth', curLineWidth);
hExcePathLossFittedCurve = plot(fittedCurve, 'k-');
set(hExcePathLossFittedCurve, 'LineWidth', curLineWidth);
hExcePathLossFittedLine = plot(fittedLine(:,1), fittedLine(:,2), 'k-.', ...
    'LineWidth', curLineWidth);
hExcePathLossFlatLine = plot([treeNumValues(1); treeNumValues(end)], ...
    [meanOfAllExcePathLosses; meanOfAllExcePathLosses], '-.', ...
    'LineWidth', curLineWidth, 'Color', ones(1,3).*0.7);
legend([hExcePathLosses, hExcePathLossMeans, hExcePathLossMedians, ...
    hExcePathLossFittedCurve, hExcePathLossFittedLine, ...
    hExcePathLossFlatLine], ...
    'Excess path losses', 'Mean', 'Median', ...
    'Fitted curve', 'Fitted line through (0,0)', 'Mean for all data', ...
    'Location','southeast');
title({'Excess Path Losses vs. Number of Trees in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excess Path Losses (dB)');
axis tight;
% Extend the x range a little bit to better show the dots on x = 0.
curAxis = axis; axis([-0.5 curAxis(2:end)]);
grid minor; transparentizeCurLegends;

saveas(hExcePathLossOverTreeNumMore, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNumMore.fig'));
saveas(hExcePathLossOverTreeNumMore, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNumMore.png'));

% What if we shift the free space model with a constant offset whenever
% there are trees present?
constShiftedFreeSpacePathLosses = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    curBoolsShiftNotNeeded = numsOfTreesInFirstFresnel{idxS}==0;
    constShiftedFreeSpacePathLosses{idxS} = freeSpacePathLosses{idxS} ...
        + meanOfAllExcePathLosses;
    constShiftedFreeSpacePathLosses{idxS}(curBoolsShiftNotNeeded) ...
        = freeSpacePathLosses{idxS}(curBoolsShiftNotNeeded);
end
allConstShiftedFreeSpacePathLosses ...
    = vertcat(constShiftedFreeSpacePathLosses{:});

% Measured path losses and free-space path losses over linear distance.
hMeasAndFreeSpaceLossesLogDist = figure;
% Free-space path losses.
hFreeSpace = semilogx(sortedLosDist, ...
    allFreeSpacePathLosses(indicesForSortedDist), '-k', 'LineWidth', 1);
hold on;
% Constantly shifted free-space path losses.
hConstShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allConstShiftedFreeSpacePathLosses(indicesForSortedDist), '*g', ...
    'MarkerSize', 4);
% Measured path losses.
hMeas = semilogx(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1), 'ob', ...
    'MarkerSize', 3);
% Shifted free-space path losses.
hShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allShiftedFreeSpacePathLosses(indicesForSortedDist), '.r', ...
    'MarkerSize', 8);
legend([hMeas, hFreeSpace, hShiftedFreeSpace, hConstShiftedFreeSpace], ...
    'Measurements', 'Free Space Path Loss', ...
    'Shifted Free Space Path Loss', ...
    'Constantly shifted Free Space Path Loss', ...
    'Location', 'southeast');
axis tight; transparentizeCurLegends;
rmseConstShiftVsMeas = sqrt(sum((allConstShiftedFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allConstShiftedFreeSpacePathLosses));
title({['Foilage Analysis Results: RMSE = ', ...
    num2str(rmseShiftVsMeas, '%.4f'), ' dB'], ...
    ['(Ref RMSE with a constant offset ', ...
    num2str(meanOfAllExcePathLosses, '%.4f'), ' dB for trees: ', ...
    num2str(rmseConstShiftVsMeas, '%.4f'), ' dB)']});
xlabel('RX Distance (m)'); ylabel('Free-Space Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLossesLogDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist.fig'));
saveas(hMeasAndFreeSpaceLossesLogDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist.png'));

% Excess path loss (z) vs tree number (y) vs distance (x).
hExcePathLossVsTreeNumVsDist = figure; hold on;
% Constantly shifted free-space path losses.
plot3k([sortedLosDist allTreeNumsInFirstFresnel(indicesForSortedDist) ...
    allExceLossRefFreeSpace(indicesForSortedDist)]);
xlabel('RX Distance (m)'); ylabel('Tree Number'); 
zlabel('Free-Space Path Loss (dB)'); grid on;

saveas(hExcePathLossVsTreeNumVsDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'ExcePathLossVsTreeNumVsDist.fig'));
saveas(hExcePathLossVsTreeNumVsDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'ExcePathLossVsTreeNumVsDist.png'));

disp('    Done!')
% EOF