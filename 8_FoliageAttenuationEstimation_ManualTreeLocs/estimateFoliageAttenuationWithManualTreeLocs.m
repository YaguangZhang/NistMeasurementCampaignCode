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

tx3D = [xTx, yTx, TX_ALT+TX_HEIGHT_M];

if exist(pathToSaveFoliageAttenResults, 'file')
    disp('        Found history results... ')
    disp('        Loading... ')
    load(pathToSaveFoliageAttenResults)
else
    [numsOfTreesInFirstFresnel, treesUtmXYHsInZone, weightedTrunkNums] ...
        = deal(cell(numOfSeries, 1));
    weightedNumsOfTreesInFirstFresnel = cell(numOfSeries, 1);
    
    for idxS = 1:numOfSeries
        disp(['        Series ', ...
            num2str(idxS), '/' num2str(numOfSeries), '...']);
        [curNumSamps, ~] = size(pathLossUtmXYHs{idxS});
        
        [numsOfTreesInFirstFresnel{idxS}, weightedTrunkNums{idxS}] ...
            = deal(nan(curNumSamps, 1));
        treesUtmXYHsInZone{idxS} = cell(curNumSamps,1);
        
        for idxSamp = 1:curNumSamps
            rx3D = pathLossUtmXYHs{idxS}(idxSamp, :);
            % Considering the RX height over the ground.
            rx3D(:, 3) = rx3D(:, 3)+RX_HEIGHT_M;
            [numsOfTreesInFirstFresnel{idxS}(idxSamp), ...
                treesUtmXYHsInZone{idxS}{idxSamp}] ...
                = countNumOfTreesInFirstFresnelZone(tx3D, rx3D, ...
                treeUtmXYHs, F_C_IN_GHZ);
            
            % Compute the weighted number of trunks in the first Fresnel
            % zone.
            [numOfTreesUtmXYHsInZone, ~] ...
                = size(treesUtmXYHsInZone{idxS}{idxSamp});
            assert( ...
                numOfTreesUtmXYHsInZone ...
                ==numsOfTreesInFirstFresnel{idxS}(idxSamp), ...
                'The number of trees in the first Fresnel zone does not match their locations!');
            
            weightedTrunkNums{idxS}(idxSamp) = computeWeightedTrunkNum( ...
                treesUtmXYHsInZone{idxS}{idxSamp}, tx3D, rx3D);
        end
    end
    
    save(pathToSaveFoliageAttenResults, ...
        'numsOfTreesInFirstFresnel', ...
        'weightedTrunkNums');
end

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

save(pathToSaveFoliageAttenResults, ...
    'freeSpacePathLosses', 'txRxLosDists', ...
    'exceLossRefFreeSpace', '-append');

disp('    Done!')

%% Plot Figures for the Tree Numbers

disp(' ')
disp('    Generating plots for tree numbers ...')

pathToSaveTreeNumPlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'TreeNums');

% Tree numbers on a 2D map.
hTreeNumOnMap = figure;
hold on;
allNumsOfTreesInFirstFresnel = vertcat(numsOfTreesInFirstFresnel{:});
allContiPathLossesWithGpsInfo = vertcat(contiPathLossesWithGpsInfo{:});
% We need (lon, lat) for plotting.
allRxGpsLocs = allContiPathLossesWithGpsInfo(:, [3 2]);
% TX.
hTx = plot3(TX_LON, TX_LAT, 0, 'kx');
% Trees.
hTrees = plot3(treeLocations(:,2), treeLocations(:,1), ...
    treeLocations(:,3).*0, 'g*');
% Tree numbers in 1st Fresnel zone for each RX location.
plot3k([allRxGpsLocs allNumsOfTreesInFirstFresnel], 'PlotType', 'stem');
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
plot(sortedLosDist, allNumsOfTreesInFirstFresnel(indicesForSortedDist));
ylabel('Number of Trees');
% Excess path losses.
yyaxis right;
plot(sortedLosDist, allExceLossRefFreeSpace(indicesForSortedDist));
title({'Number of Trees in the 1st Fresnel Zone', ...
    'with Excess Path Losses (Over Free-Space Losses)'});
xlabel('Distance to TX (m)'); ylabel('Excess Path Loss (dB)'); grid on;

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
xlabel('Distance to TX (m)'); ylabel('Free-Space Path Loss (dB)'); grid on;

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
xlabel('Distance to TX (m)'); ylabel('Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLossesSameY, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.fig'));
saveas(hMeasAndFreeSpaceLossesSameY, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'MeasAndFreeSpaceLosses.png'));

% Excess path loss over number of trees in the 1st Fresnel zone.
hExcePathLossOverTreeNum = figure;
plot(allNumsOfTreesInFirstFresnel, allExceLossRefFreeSpace, '.');
title({'Excess Path Losses vs. Number of Trunks in the 1st Fresnel Zone'});
xlabel('Number of Trees'); ylabel('Excess Path Losses (dB)');
axis tight; grid on;

saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.fig'));
saveas(hExcePathLossOverTreeNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNum.png'));

% Excess path loss over log number of trees in the 1st Fresnel zone.
hExcePathLossOverLogTreeNum = figure;
semilogx(allNumsOfTreesInFirstFresnel, allExceLossRefFreeSpace, '.');
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

% The excess path loss vs weighted trunk number.
allWeightedTrunkNums = vertcat(weightedTrunkNums{:});

hExcePathLossOverWeightedTrunkNum = figure; hold on;
scatter(allWeightedTrunkNums, allExceLossRefFreeSpace, '.');
xlabel('Weighted trunk number'); ylabel('Excess Path Losses (dB)');
axis tight; grid minor;

saveas(hExcePathLossOverWeightedTrunkNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'ExcePathLossOverWeightedTrunkNum.fig'));
saveas(hExcePathLossOverWeightedTrunkNum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'ExcePathLossOverWeightedTrunkNum.png'));

disp('    Done!')

%% Apply the Partition-Dependent Model
% Essentially we will get one fixed value for the extra path loss caused by
% each tree inside the 1st Fresnel zone.

disp(' ')
disp('    Computing the extra path loss caused by each tree in the 1st Fresnel zone...')

excePathLossPerTree ...
    = pinv(allNumsOfTreesInFirstFresnel)*allExceLossRefFreeSpace;

disp('    Making predictions accordingly ...')
predictionsConstLossPerTrunk = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);
    
    predictionsConstLossPerTrunk{idxS} = freeSpacePathLosses{idxS} ...
        + numsOfTreesInFirstFresnel{idxS}.*excePathLossPerTree;
end

save(pathToSaveFoliageAttenResults, 'excePathLossPerTree', ...
    'predictionsConstLossPerTrunk', '-append');

% What if we shift the free space model with constant offsets for each
% tree, but with the constant offsets estimated within each tree number
% group?

disp(' ')
disp('    Computing the extra path loss for each tree number group separately...')

% We have on each row: group number, number of trees, mean extra path loss
% in dB (for locations in that tree number group).
allNumsOfTreesInFirstFresnel = vertcat(numsOfTreesInFirstFresnel{:});
numOfTreeNumberGroups = max(allNumsOfTreesInFirstFresnel)...
    -min(allNumsOfTreesInFirstFresnel)+1;
excePathLossGroupWise = nan(numOfTreeNumberGroups, 3);

% Compute the average excessive path loss in each tree number group.
treeNumValues = min(allNumsOfTreesInFirstFresnel) ...
    :max(allNumsOfTreesInFirstFresnel);
excePathLossMean = arrayfun(@(n) ...
    mean(allExceLossRefFreeSpace(allNumsOfTreesInFirstFresnel==n)), ...
    treeNumValues);
assert(length(excePathLossMean)==numOfTreeNumberGroups, ...
    'There should be one excessive path loss value computed for each tree number group!')

for idxTreeNumG = 1:numOfTreeNumberGroups
    curNumTrees = min(allNumsOfTreesInFirstFresnel)+idxTreeNumG-1;
    curExcePathLoss = excePathLossMean(idxTreeNumG);
    excePathLossGroupWise(idxTreeNumG,:) ...
        = [idxTreeNumG, curNumTrees, curExcePathLoss];
end

disp('    Again, making predictions accordingly ...')
groupWiseShiftedFreeSpacePathLosses = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    disp(['        Series ', ...
        num2str(idxS), '/' num2str(numOfSeries), '...']);
    
    groupWiseExcePathLoss = arrayfun(@(nTrees) ...
        excePathLossGroupWise( ...
        excePathLossGroupWise(:,2)==nTrees, 3), ...
        numsOfTreesInFirstFresnel{idxS});
    groupWiseShiftedFreeSpacePathLosses{idxS} = freeSpacePathLosses{idxS} ...
        + groupWiseExcePathLoss;
end

save(pathToSaveFoliageAttenResults, 'excePathLossGroupWise', ...
    'groupWiseShiftedFreeSpacePathLosses', '-append');

%% Plot the Path Loss Results

allPredictionsConstLossPerTrunk = vertcat(predictionsConstLossPerTrunk{:});
allGroupWiseShiftedFreeSpacePathLosses ...
    = vertcat(groupWiseShiftedFreeSpacePathLosses{:});

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
    allPredictionsConstLossPerTrunk(indicesForSortedDist), '-');
legend([hMeas, hFreeSpace, hShiftedFreeSpace], ...
    'Measurements', 'Free Space Path Loss', ...
    'Shifted Free Space Path Loss', 'Location', 'southeast');
axis tight; transparentizeCurLegends;

% RMSE
rmseShiftVsMeas = sqrt(sum((allPredictionsConstLossPerTrunk ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allPredictionsConstLossPerTrunk));
rmseGroupWiseShiftVsMeas = sqrt(...
    sum((allGroupWiseShiftedFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allGroupWiseShiftedFreeSpacePathLosses));

title(['Foilage Analysis Results: RMSE = ', ...
    num2str(rmseShiftVsMeas, '%.4f'), ' dB']);
xlabel('Distance to TX (m)'); ylabel('Free-Space Path Loss (dB)'); grid on;

saveas(hMeasAndFreeSpaceLossesLinearDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLinearDist.fig'));
saveas(hMeasAndFreeSpaceLossesLinearDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLinearDist.png'));

% Excess path loss over number of trees in the 1st Fresnel zone with some
% analysis information overlaid.
hExcePathLossOverTreeNumMore = figure; hold on;
numPtsToShow = length(allNumsOfTreesInFirstFresnel);
hExcePathLosses ...
    = scatter(allNumsOfTreesInFirstFresnel, allExceLossRefFreeSpace, ...
    ones(numPtsToShow, 1) .* 8, ...
    'MarkerFaceColor','b','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% Also show the median, mean and fitted lines.
excePathLossMedian = arrayfun(@(n) ...
    median(allExceLossRefFreeSpace(allNumsOfTreesInFirstFresnel==n)), ...
    treeNumValues);
fittedCurve = fit(allNumsOfTreesInFirstFresnel, allExceLossRefFreeSpace, ...
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

% What if we shift the free space model with a constant offset?
constShiftedFreeSpacePathLosses = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    constShiftedFreeSpacePathLosses{idxS} = freeSpacePathLosses{idxS} ...
        + meanOfAllExcePathLosses;
end
allConstShiftedFreeSpacePathLosses ...
    = vertcat(constShiftedFreeSpacePathLosses{:});

% What if we shift the free space model with a constant offset whenever
% there are trees present? Similar to allConstShiftedFreeSpacePathLosses,
% just this time we are interested in the average offset for only locations
% with tree blockages.
excePathLossesWithTrees = arrayfun(@(idxS) ...
    exceLossRefFreeSpace{idxS}(numsOfTreesInFirstFresnel{idxS}>0), ...
    1:numOfSeries, 'UniformOutput', false);
meanOfExcePathLossesWithTrees = mean(vertcat(excePathLossesWithTrees{:}));
constShiftedForTreesFreeSpacePathLosses = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    curBoolsShiftNotNeeded = numsOfTreesInFirstFresnel{idxS}==0;
    constShiftedForTreesFreeSpacePathLosses{idxS} = freeSpacePathLosses{idxS} ...
        + meanOfExcePathLossesWithTrees;
    constShiftedForTreesFreeSpacePathLosses{idxS}(curBoolsShiftNotNeeded) ...
        = freeSpacePathLosses{idxS}(curBoolsShiftNotNeeded);
end
allConstShiftedForTreesFreeSpacePathLosses ...
    = vertcat(constShiftedForTreesFreeSpacePathLosses{:});

% Plot the models with the measurement results.
hMeasAndFreeSpaceLossesLogDist = figure;
% Widen the figure.
with2HeightRatio = 2;
curPos = get(hMeasAndFreeSpaceLossesLogDist, 'Pos');
curPos(3) = curPos(4) *with2HeightRatio;
set(hMeasAndFreeSpaceLossesLogDist, 'Pos', curPos);
% Free-space path losses.
hFreeSpace = semilogx(sortedLosDist, ...
    allFreeSpacePathLosses(indicesForSortedDist), '-', ...
    'Color', ones(1,3).*0.7, 'LineWidth', 1);
hold on;
% Constantly shifted free-space path losses.
hConstShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allConstShiftedFreeSpacePathLosses(indicesForSortedDist), ...
    '*m', 'MarkerSize', 4);
% Constantly shifted free-space path losses, only for locations blocked by
% trees.
hConstShiftedForTreeBlockageFreeSpace ...
    = semilogx(sortedLosDist, ...
    allConstShiftedForTreesFreeSpacePathLosses(indicesForSortedDist), ...
    '^g', 'MarkerSize', 4);
% Measured path losses.
hMeas = semilogx(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1), 'ob', ...
    'MarkerSize', 3);
% Shifted free-space path losses, with constant extra path loss per tree.
hShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allPredictionsConstLossPerTrunk(indicesForSortedDist), '.r', ...
    'MarkerSize', 6);
% Shifted free-space path losses, with constant extra path loss per tree
% for locations in each tree number group.
hGroupWiseShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allGroupWiseShiftedFreeSpacePathLosses(indicesForSortedDist), '.', ...
    'Color', ones(1,3).*0.5, 'MarkerSize', 6);
hCurLegend = legend([hMeas, hFreeSpace, ...
    hShiftedFreeSpace, hGroupWiseShiftedFreeSpace, ...
    hConstShiftedForTreeBlockageFreeSpace, hConstShiftedFreeSpace], ...
    'Measurements', ...
    'Free-space path loss (FSPL)', ...
    'FSPL with a constant loss per tree', ...
    'FSPL with an extra group-wise loss', ...
    'FSPL with an extra NLoS loss', ...
    'FSPL with an extra loss', ...
    'Location', 'northwest');
% Manually adjust the legend postion to avoid blocking data points.
% set(hCurLegend, 'Pos', [0.1405, 0.6367, 0.4278, 0.2726]);
axis tight; transparentizeCurLegends;
rmseConstShiftVsMeas = sqrt(sum((allConstShiftedFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allConstShiftedFreeSpacePathLosses));
% title({['Foilage Analysis Results: RMSE = ', ...
%     num2str(rmseShiftVsMeas, '%.4f'), ' dB'], ...
%      ['(Ref RMSE with a constant offset ', ...
%     num2str(meanOfAllExcePathLosses, '%.4f'), ' dB for trees: ', ...
%      num2str(rmseConstShiftVsMeas, '%.4f'), ' dB)']});
xlabel('Distance to TX (m)'); ylabel('Path Loss (dB)'); grid on;
% Manually set the x tick labels.
xticks([10, 20, 30, 50, 100, 200]);
xticklabels(arrayfun(@(x) {num2str(x)}, xticks'))
saveas(hMeasAndFreeSpaceLossesLogDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist.fig'));
saveas(hMeasAndFreeSpaceLossesLogDist, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist.png'));

% Export an .eps copy for papers.
pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', '1_EpsFigs');
saveEpsFigForPaper(hMeasAndFreeSpaceLossesLogDist, ...
    fullfile(pathToSavePaperFigs, ...
    ['4_', 'MeasAndFreeSpaceLossesAndShiftedOverLogDist.eps']));

% Generate a simpler version of this model comparison plot. Plot the models
% with the measurement results.
hMeasAndFreeSpaceLossesLogDistSimp = figure;
% Free-space path losses.
hFreeSpace = semilogx(sortedLosDist, ...
    allFreeSpacePathLosses(indicesForSortedDist), '-', ...
    'Color', ones(1,3).*0.7, 'LineWidth', 1);
hold on;
% Measured path losses.
hMeas = semilogx(sortedLosDist, ...
    allContiPathLossesWithGpsInfo(indicesForSortedDist,1), 'ob', ...
    'MarkerSize', 3);
% Constantly shifted free-space path losses.
hConstShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allConstShiftedFreeSpacePathLosses(indicesForSortedDist), ...
    '*g', 'MarkerSize', 4);
% Shifted free-space path losses, with constant extra path loss per tree.
hShiftedFreeSpace ...
    = semilogx(sortedLosDist, ...
    allPredictionsConstLossPerTrunk(indicesForSortedDist), '.r', ...
    'MarkerSize', 8);
hCurLegend = legend([hMeas, hFreeSpace, ...
    hShiftedFreeSpace, ...
    hConstShiftedFreeSpace], ...
    'Measurements', ...
    'Free-space path loss (FSPL)', ...
    'FSPL with a constant loss per tree', ...
    'FSPL with an extra loss', ...
    'Location', 'northwest');
axis tight; transparentizeCurLegends;
xlabel('Distance to TX (m)'); ylabel('Path Loss (dB)'); grid on;
% Manually set the x tick labels.
xticks([10, 20, 30, 50, 100, 200]);
xticklabels(arrayfun(@(x) {num2str(x)}, xticks'));
saveas(hMeasAndFreeSpaceLossesLogDistSimp, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist_Simplified.fig'));
saveas(hMeasAndFreeSpaceLossesLogDistSimp, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist_Simplified.png'));

% Export an .eps copy for papers.
pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', '1_EpsFigs');
saveEpsFigForPaper(hMeasAndFreeSpaceLossesLogDistSimp, ...
    fullfile(pathToSavePaperFigs, ...
    ['4_1_', ...
    'MeasAndFreeSpaceLossesAndShiftedOverLogDist_Simplified.eps']));

% Save RMSE results to a txt file.
rmseFsplVsMeas = sqrt(sum((allFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allFreeSpacePathLosses));
rmseConstShiftForTreesVsMeas = sqrt(sum(( ...
    allConstShiftedForTreesFreeSpacePathLosses ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allConstShiftedForTreesFreeSpacePathLosses));

disp('    Done!')

%% Plots for Publication

disp(' ')
disp('    Generating plots for publication ...')

% Modify hExcePathLossOverTreeNumMore to generate a figure for publication.
hExcePathLossOverTreeNumModels = figure;
curFigPos = get(gcf, 'Position');
set(gcf, 'Position', [curFigPos(1) 0 ...
    ones(1,2).*max(curFigPos(3:4))+[0, -60]]);

% hSupTitle = suptitle( ...
%     {'Excess Path Losses vs. Number of Trees in the 1st Fresnel Zone'});
% set(hSupTitle, 'FontSize', 12);

% Set this to arrange the subfigures.
flagModelCompOnTop = true;
% Model comparison.
if flagModelCompOnTop
    subplot('Position',[.1 .225 .85 .675]);
else
    subplot('Position',[.1 .1 .85 .675]);
end
hold on;

% One can use
%   boolsTreeNumValuesToShow = treeNumValues>0;
% and
%   boolsTreeNumsInFirstFresnelToShow = allTreeNumsInFirstFresnel>0;
% to show only data for positive tree numbers.

boolsTreeNsInFirstFresToShow = allNumsOfTreesInFirstFresnel>0;
boolsTreeNumsToShow = treeNumValues>0;
treeNumValuesToShow = treeNumValues(boolsTreeNumsToShow);
% Add an error bar to better show the distribution of the path loss
% measurement results.
excePathLossStD = arrayfun(@(n) ...
    std(allExceLossRefFreeSpace(allNumsOfTreesInFirstFresnel==n)), ...
    treeNumValuesToShow);

hErrorBar = errorbar(treeNumValues(boolsTreeNumsToShow), ...
    excePathLossMean(boolsTreeNumsToShow), ...
    excePathLossStD, 'LineStyle', 'none', ...
    'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);
% Measurement results.
hExcePathLosses ...
    = scatter(allNumsOfTreesInFirstFresnel(boolsTreeNsInFirstFresToShow), ...
    allExceLossRefFreeSpace(boolsTreeNsInFirstFresToShow), ...
    ones(sum(boolsTreeNsInFirstFresToShow), 1) .* 8, ...
    'MarkerFaceColor','b','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% Also show the median.
excePathLossMedian = arrayfun(@(n) ...
    median(allExceLossRefFreeSpace(allNumsOfTreesInFirstFresnel==n)), ...
    treeNumValuesToShow);
g = fittype('a+x^b');
modExpPoly = fit(allNumsOfTreesInFirstFresnel(boolsTreeNsInFirstFresToShow), ...
    allExceLossRefFreeSpace(boolsTreeNsInFirstFresToShow), g, ...
    'StartPoint', rand(1,2));

% The fitted line will be stored as [x(1), y(1); x(2) y(2)].
fittedLine = [[0; treeNumValues(end)], ...
    [0; treeNumValues(end)].*excePathLossPerTree];
% Also plot a flat line (a constant offset with slope 0).
meanOfAllExcePathLosses ...
    = mean(allExceLossRefFreeSpace(boolsTreeNsInFirstFresToShow));
% Plots.
curLineWidth = 1.5;

hExcePathLossMedians = plot(treeNumValues(boolsTreeNumsToShow), ...
    excePathLossMedian, 'xb--', ...
    'LineWidth', curLineWidth);

xx = linspace(0,max(treeNumValues)+1,100);
hExcePathLossFittedCurve = plot(xx, modExpPoly(xx), 'k-');
set(hExcePathLossFittedCurve, 'LineWidth', curLineWidth);

hExcePathLossFlatLine = plot(...
    [treeNumValues(1); treeNumValues(end)+1], ...
    [meanOfAllExcePathLosses; meanOfAllExcePathLosses], '-.', ...
    'LineWidth', curLineWidth, 'Color', ones(1,3).*0.7);

legend([hExcePathLosses, hErrorBar, hExcePathLossMedians, ...
    hExcePathLossFittedCurve, ...
    hExcePathLossFlatLine], ...
    'Excess path losses', 'One standard deviation', 'Median', ...
    'Exponential polynomial (a+x^b)', 'Mean for all data', ...
    'Location','northeast');
if flagModelCompOnTop
    xticklabels({});
else
    xlabel('Number of Trees');
end
ylabel('Excess Path Losses (dB)');
axis tight;
% Extend the x range a little bit to better show the dots on x = 0.
curAxis = axis; deltaY = (curAxis(4)-curAxis(3)).*0.05;
axis([curAxis(1) curAxis(2) ...
    curAxis(3)-deltaY curAxis(4)+3.*deltaY]);
curAxis = axis;
grid on; grid minor; transparentizeCurLegends;

% Bar chart for number of measurement results.
if flagModelCompOnTop
    hBarChart = subplot('Position',[.1 .1 .85 .125]);
else
    hBarChart = subplot('Position',[.1 .775 .85 .125]);
end
hold on;
ax = ancestor(hBarChart, 'axes');
ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;

excePathLossFreq = arrayfun(@(n) ...
    sum(allNumsOfTreesInFirstFresnel==n), ...
    treeNumValuesToShow);
bar(treeNumValuesToShow, excePathLossFreq);
curAxisBar = axis; deltaY = curAxisBar(4) - curAxisBar(3);
axis([curAxis(1:2) curAxisBar(3) curAxisBar(4)+deltaY/10]);
grid on; grid minor;
if flagModelCompOnTop
    xlabel('Number of Trunks in the 1st Fresnel Zone');
else
    xticklabels({});
end

% ylabel({'Number of','Samples'});
hBarChartLeg = legend('Number of path loss results');
set(hBarChartLeg, 'FontSize', 9);

% Save the figure.
saveas(hExcePathLossOverTreeNumModels, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNumModels.fig'));
saveas(hExcePathLossOverTreeNumModels, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ExceLossOverTreeNumModels.png'));

% Export an .eps copy for papers.
pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', '1_EpsFigs');
saveEpsFigForPaper(hExcePathLossOverTreeNumModels, ...
    fullfile(pathToSavePaperFigs, ...
    ['3_2_', 'ExceLossOverTreeNumModels.eps']));
% Also save a .png copy there because the overlapping effect of transparent
% markers is lost in the .eps figure.
saveas(hExcePathLossOverTreeNumModels, ...
    fullfile(pathToSavePaperFigs, ...
    ['3_2_', 'ExceLossOverTreeNumModels.png']));

% RMSE
shiftedFsplModExp = cell(numOfSeries, 1);
for idxS = 1:numOfSeries
    shiftedFsplModExp{idxS} = freeSpacePathLosses{idxS} ...
        + modExpPoly(numsOfTreesInFirstFresnel{idxS});
end
allShiftedFsplModExp = vertcat(shiftedFsplModExp{:});
rmseModExpDecVsMeas = sqrt(sum((allShiftedFsplModExp ...
    - allContiPathLossesWithGpsInfo(:,1)).^2)...
    /length(allShiftedFsplModExp));

%% Save MSEs to a file.
dbLossFormatter = '%3.2f';
fIdRmse = fopen( ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'RmseComparisons.txt'), ...
    'w');
fprintf(fIdRmse,'Model, Extra Loss(dB), RMSE(dB)\n');
% FSPL + group-wise const extra loss per tree
fprintf(fIdRmse, ['FSPL with an extra group-wise loss, ', ...
    'different for each tree number group, ', dbLossFormatter, '\r'], ...
    rmseGroupWiseShiftVsMeas);
% FSPL + exponential polynomial (a+x^b)
fprintf(fIdRmse, ['FSPL with an exponential polynomial (a+x^b) decay, ', ...
    'different for each tree number group, ', dbLossFormatter, '\r'], ...
    rmseModExpDecVsMeas);
% FSPL + const extra loss
fprintf(fIdRmse, ['FSPL with an extra loss, ', ...
    dbLossFormatter, ', ', dbLossFormatter, '\r'], ...
    meanOfAllExcePathLosses, rmseConstShiftVsMeas);
% FSPL + const extra loss for foliage blockages
fprintf(fIdRmse, ['FSPL with an extra NLoS loss, ', ...
    dbLossFormatter, ', ', dbLossFormatter, '\r'], ...
    meanOfExcePathLossesWithTrees, rmseConstShiftForTreesVsMeas);
% FSPL + const extra loss per tree
fprintf(fIdRmse, ['FSPL with a constant loss for each tree, ', ...
    dbLossFormatter, ', ', dbLossFormatter, '\r'], ...
    excePathLossPerTree, rmseShiftVsMeas);
% FSPL
fprintf(fIdRmse, ['FSPL, 0, ', ...
    dbLossFormatter, '\r'], ...
    rmseFsplVsMeas);
fclose(fIdRmse);

disp('    Done!')

% EOF