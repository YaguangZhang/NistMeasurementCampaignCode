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

% Reuse results from plotInfo.m, calibrateRx.m, fetchAntennaPattern.m, and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_PATH_LOSSES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');

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
disp('      - contiPathLossesWithGpsInfo.m')

assert(exist(ABS_PATH_TO_PATH_LOSSES, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/3_PathLossComputation/evalPathLossesForContiTracks.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'contiPathLossesWithGpsInfo', 'contiOutFilesRelPathsUnderDataFolder'
% and 'contiOutFileIndicesReflection'.
load(ABS_PATH_TO_PATH_LOSSES);

disp('    Done!')

%% Load the Tree Info

disp(' ')
disp('    Loading tree location information...')

treeLocations = loadGpsMarkers(ABS_PATH_TO_TREE_LOCS, 'Marker*', ...
    pathToSaveTreeLocs);
% Keep trying to fetch alts from Google until all needed info is available.
numRetrial = 0;
while any(isnan(treeLocations(:,3)))
    numRetrial = numRetrial+1;
    disp(['Not all alts retrieved. Retry again ... (#', ...
        num2str(numRetrial), ')']);
    treeLocations = loadGpsMarkers(ABS_PATH_TO_TREE_LOCS, 'Marker*', ...
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

% We are not yet able to deal with points from different UTM zones.
allRelevantZones ...
    = num2cell([vertcat(pathLossUtmZones{:});treeUtmZones], 2);
firstZone = allRelevantZones(1,:);
assert(all(strcmp(allRelevantZones, firstZone)), ...
    'We are not yet able to deal with points from different UTM zones!');

save(pathToSaveUtmInfo, ...
    'pathLossUtmXYHs', 'pathLossUtmZones', ...
    'treeUtmXYHs', 'treeUtmZones');

disp('    Done!')

%% Plot an Overview in UTM

pathToSaveTreeLocs = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'OverviewInUtm');
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

saveas(hOverviewInUtm, [pathToSaveTreeLocs, '.fig']);
saveas(hOverviewInUtm, [pathToSaveTreeLocs, '.png']);

%% Count the Number of Trees in the Fresii for Each Path Loss Measurement



% EOF