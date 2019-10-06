% INSPECTOUTLAYERS Inspect outlayers in the path loss measurements to
% investigate the cause.
%
% Yaguang Zhang, Purdue, 07/29/2019

clc; close all;

% Locate the current working directory.
cd(fileparts(mfilename('fullpath')));
[~,folderNameToSaveResults,~] = fileparts(pwd);
cd('..'); addpath('lib');
curFileName = mfilename;
fileNameHintRuler = hintScriptName(curFileName);

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
setPath;

%% Before Processing the Data

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', folderNameToSaveResults);

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% We will also generate plots to show the half-power beamwidth (HPBW) and
% first-null beamwidth (FNBW) of the TX.
azimuthHpbwInDeg = 10.08;
azimuthFnbwInDeg = 29.716+(360-331.961);

%% Load Measurement Data

ABS_PATH_TO_TX_INFO_LOGS_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');
ABS_PATH_TO_PATH_LOSSES_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti', ...
    'contiPathLossesWithGpsInfo.mat');
ABS_PATH_TO_MODEL_PREDICTIONS_FILE ...
    = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', ...
    'FoliageDepthBasedModelsComparison', ...
    'predictionsFromSelectedModels.mat');
try
    % Get records of the TxInfo.txt files (among other contant parameters
    % for the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and
    % TX_POWER_DBM):
    % 	'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
    load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
    % Get 'contiPathLossesWithGpsInfo',
    % 'contiOutFilesRelPathsUnderDataFolder' and
    % 'contiOutFileIndicesReflection'.
    load(ABS_PATH_TO_PATH_LOSSES_FILE);
catch
    error('Unable to load the measurement data!');
end

% Extract the data for readability.
txLat = TX_LAT;
txLon = TX_LON;

allContiPathLossesWithGpsInfo = vertcat(contiPathLossesWithGpsInfo{:});
rxLats = allContiPathLossesWithGpsInfo(:, 2);
rxLons = allContiPathLossesWithGpsInfo(:, 3);

% Any prediction error with absolute value beyond this value will be
% considered as an outlayer.
MAX_ALLOWED_ERROR_IN_DB = 3;

%% Load Model Predictions
% We need to compare the predictions with the measurements (rxPathlosses)
% and generate plots for each model listed below:
%   - Free space path loss (FSPL) model => rxFSPLs
%    - Attenuation factor (AF) model
%   - Weissberger's modified exponential decay (WMED) model
%    - ITU-R obstruction by woodland model
%   - Site-specific model A-I
%    - Site-specific model A-II
%   - Site-specific model B
%    - Site-specific model C

% Measurements.
allContiPathLossesWithGpsInfo = vertcat(contiPathLossesWithGpsInfo{:});
try
    % Get 'allFreeSpacePathLosses', 'allPredictionsConstLossPerTrunk',
    % 'allPredictedPathLossesItu', 'allPredictedPathLossesMod',
    % 'allPredictedPathLossesTwoStepConFixedB',
    % 'allPredictedPathLossesTwoStepCLPerUoFAndCLForDEV',
    % 'allPredictedPathLossesItuModForFoliageDepth', and
    % 'allPredictedPathLossesTwoStepLinearLossWrtFA'.
    load(ABS_PATH_TO_MODEL_PREDICTIONS_FILE);
catch
    error('Not able to find the model prediction file!');
end

% Rename the data for readability.
allMeas = allContiPathLossesWithGpsInfo(:, 1);
measLegend = 'Measurements';

modelsToInspect = {'FSPL', allFreeSpacePathLosses; ...
    'AF', allPredictionsConstLossPerTrunk; ...
    'ITU', allPredictedPathLossesItu; ...
    'WMED', allPredictedPathLossesMod; ...
    'A-I', allPredictedPathLossesTwoStepConFixedB; ...
    'A-II', allPredictedPathLossesTwoStepCLPerUoFAndCLForDEV; ...
    'B', allPredictedPathLossesItuModForFoliageDepth; ...
    'C', allPredictedPathLossesTwoStepLinearLossWrtFA};
modelLegends = {'Free space path loss', ...
    'Constant loss per tree', ...
    'ITU obstruction by woodland', ...
    "Weissberger's model", ...
    'Site-specific model A-I', 'Site-specific model A-II', ...
    'Site-specific model B', 'Site-specific model C'};
numModelsToInspect = length(modelLegends);

%% Figs

FLAG_GEN_FIGS_SILENTLY = true;

for idxModel = 1:numModelsToInspect
    curModelName = modelsToInspect{idxModel,1};
    curDirToSaveFigs = fullfile(pathToSaveResults, ...
        ['errorsOnMap_', curModelName]);
    
    curModelPredictions = modelsToInspect{idxModel,2};
    curModelLegend = modelLegends{idxModel};
    
    curErrors = curModelPredictions-allMeas;
    
    rxLocMarker = 'og';
    rxLocMarkerSize = 4;
    rxLocMarkerFaceColor = 'none';
    rxLocMarkerLineWidth = 0.15;
    
    numOfTicklabels = 11;
    
    % Fig 1: absolute errors on map
    hFigErrorsOnMap = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    hRx = plot3(rxLons, rxLats, zeros(size(rxLats)), rxLocMarker, ...
        'MarkerFaceColor', rxLocMarkerFaceColor, ...
        'MarkerSize', rxLocMarkerSize, 'LineWidth', rxLocMarkerLineWidth);
    % Make sure the color map is intuitive and works in grey scale.
    colormap hot;
    % Shift the lowest point to +1dB.
    curAbsErrors = abs(curErrors);
    [~, ~, hPlot3kColorBar] ...
        = plot3k([rxLons, rxLats, curAbsErrors], ...
        'Marker', {'.', 12}, ...
        'ColorRange', [max(curAbsErrors), min(curAbsErrors)], ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Absolute Errors on Map for ', curModelName], ...
        'Longitude', 'Latitude', ...
        '', 'Error (dB)'});
    xticks([]); yticks([]);
    
    % Manually adjust the visible area.
    figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
        39.9893839683981, 39.9915745444857];
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = findall( allchild(hFigErrorsOnMap), 'type', 'colorbar');
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    
    hold off; view(2);
    legend([hTx, hRx], {'TX', 'RX Locs'}, 'Location','southeast');
    
    % Repeat for a better output result: Manually adjust the visible area.
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    saveas(hFigErrorsOnMap, [curDirToSaveFigs, '_absErr.jpg']);
    
    % Fig 2: positive outlayer errors on map
    hFigOutlayerErrorsOnMap = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY);
    hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    hRx = plot3(rxLons, rxLats, zeros(size(rxLats)), rxLocMarker, ...
        'MarkerFaceColor', rxLocMarkerFaceColor, ...
        'MarkerSize', rxLocMarkerSize, 'LineWidth', rxLocMarkerLineWidth);
    % Make sure the color map is intuitive and works in grey scale.
    colormap hot;
    boolsOutlayers = curErrors>MAX_ALLOWED_ERROR_IN_DB;
    
    outLayerErrors = curErrors(boolsOutlayers);
    outLayerRxLons = rxLons(boolsOutlayers);
    outLayerRxLats = rxLats(boolsOutlayers);
    
    plot3k([outLayerRxLons, outLayerRxLats, ...
        outLayerErrors], ...
        'Marker', {'.', 12}, ...
        'ColorRange', [max(outLayerErrors), min(outLayerErrors)], ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Positive Outlayer Errors (>', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), ' dB) on Map for ', ...
        curModelName], ...
        'Longitude', 'Latitude', ...
        '', 'Error (dB)'});
    xticks([]); yticks([]);
    
    % Manually adjust the visible area.
    figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
        39.9893839683981, 39.9915745444857];
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = findall( allchild(hFigOutlayerErrorsOnMap), 'type', 'colorbar');
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    
    hold off; view(2);
    legend([hTx, hRx], {'TX', 'RX Locs'}, 'Location','southeast');
    
    % Repeat for a better output result: Manually adjust the visible area.
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    saveas(hFigOutlayerErrorsOnMap, ...
        [curDirToSaveFigs, '_positiveOutlayerErr_', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), 'dB.jpg']);
    
    % Fig 3: negative outlayer errors on map
    hFigOutlayerErrorsOnMap = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY);
    hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    hRx = plot3(rxLons, rxLats, zeros(size(rxLats)), rxLocMarker, ...
        'MarkerFaceColor', rxLocMarkerFaceColor, ...
        'MarkerSize', rxLocMarkerSize, 'LineWidth', rxLocMarkerLineWidth);
    % Make sure the color map is intuitive and works in grey scale.
    colormap hot;
    boolsOutlayers = curErrors<-MAX_ALLOWED_ERROR_IN_DB;
    
    outLayerErrors = curErrors(boolsOutlayers);
    outLayerRxLons = rxLons(boolsOutlayers);
    outLayerRxLats = rxLats(boolsOutlayers);
    
    % Shift the lowest point to +1dB.
    errorShift = 1-min(outLayerErrors);
    [~, ~, hPlot3kColorBar] ...
        = plot3k([outLayerRxLons, outLayerRxLats, ...
        outLayerErrors+errorShift], ...
        'Marker', {'.', 12}, ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Negative Outlayer Errors (<-', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), ' dB) on Map for ', ...
        curModelName], 'Longitude', 'Latitude', ...
        '', 'Error (dB)'});
    xticks([]); yticks([]);
    
    % Manually adjust the visible area.
    figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
        39.9893839683981, 39.9915745444857];
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = hPlot3kColorBar;
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    
    % Update the colorbar tick labels.
    hCb.TickLabels = arrayfun(@(n) num2str(n-errorShift, '%.1f'), ...
        cellfun(@(l) str2num(l), hCb.TickLabels), 'UniformOutput', false);
    
    hold off; view(2);
    legend([hTx, hRx], {'TX', 'RX Locs'}, 'Location','southeast');
    
    % Repeat for a better output result: Manually adjust the visible area.
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    saveas(hFigOutlayerErrorsOnMap, ...
        [curDirToSaveFigs, '_negativeOutlayerErr_', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), 'dB.jpg']);
    
    % Fig 4: non-outlayer errors on map
    hFigNonOutlayerErrorsOnMap ...
        = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    hRx = plot3(rxLons, rxLats, zeros(size(rxLats)), rxLocMarker, ...
        'MarkerFaceColor', rxLocMarkerFaceColor, ...
        'MarkerSize', rxLocMarkerSize, 'LineWidth', rxLocMarkerLineWidth);
    % Make sure the color map is intuitive and works in grey scale.
    colormap hot;
    boolsNonOutlayers = curErrors>=-MAX_ALLOWED_ERROR_IN_DB ...
        & curErrors<=MAX_ALLOWED_ERROR_IN_DB;
    
    nonOutLayerErrors = curErrors(boolsNonOutlayers);
    nonOutLayerRxLons = rxLons(boolsNonOutlayers);
    nonOutLayerRxLats = rxLats(boolsNonOutlayers);
    
    % Shift the lowest point to +1dB.
    errorShift = 1-min(nonOutLayerErrors);
    [~, ~, hPlot3kColorBar] ...
        = plot3k([nonOutLayerRxLons, nonOutLayerRxLats, ...
        nonOutLayerErrors+errorShift], ...
        'Marker', {'.', 12}, ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Non-Outlayer Errors (', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), ' dB) on Map for ', ...
        curModelName], 'Longitude', 'Latitude', ...
        '', 'Error (dB)'});
    xticks([]); yticks([]);
    
    % Manually adjust the visible area.
    figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
        39.9893839683981, 39.9915745444857];
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = hPlot3kColorBar;
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    
    % Update the colorbar tick labels.
    hCb.TickLabels = arrayfun(@(n) num2str(n-errorShift, '%.1f'), ...
        cellfun(@(l) str2num(l), hCb.TickLabels), 'UniformOutput', false);
    
    hold off; view(2);
    legend([hTx, hRx], {'TX', 'RX Locs'}, 'Location','southeast');
    
    % Repeat for a better output result: Manually adjust the visible area.
    axis(figAxisToSet);
    plot_google_map('MapType','satellite');
    
    saveas(hFigNonOutlayerErrorsOnMap, ...
        [curDirToSaveFigs, '_nonOutlayerErr_', ...
        num2str(MAX_ALLOWED_ERROR_IN_DB), 'dB.jpg']);
end

%% TX Antenna Main Lobe Figs
% For each track's measurement results, generate an overview plot with the
% TX antenna coverage area. We will show both the HPBW and FNBW areas as
% triangle polyshapes.

BEAM_AREA_SIDE_LENGTH_IN_M = 500;
EXPECTED_PATH_LOSS_RANGE = [195, 50]; % [max(allMeas), min(allMeas)]
numOfTracks = length(contiPathLossesWithGpsInfo);

% Overview.
allRxLonLats = allContiPathLossesWithGpsInfo(:, [3,2]);

curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
hTx = plot(TX_LON, TX_LAT, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
colormap hot;
plot3k([allRxLonLats, allMeas], ...
    'Marker', {'.', 12}, ...
    'ColorRange', EXPECTED_PATH_LOSS_RANGE, ...
    'CBLabels', numOfTicklabels, 'Labels', ...
    {'Overview for Path Loss', ...
    'Longitude', 'Latitude', ...
    '', 'Path Loss (dB)'});
xticks([]); yticks([]);
axis(figAxisToSet); view(2);
plotGoogleMapAfterPlot3k(curFig, 'satellite');

curDirToSaveFig = fullfile(pathToSaveResults, ...
    ['Overview_BasicPLs.png']);
saveas(curFig, curDirToSaveFig);

% By track.
for idxTrack = 1:numOfTracks
    curPLs = contiPathLossesWithGpsInfo{idxTrack}(:,1);
    curRxLonLats = contiPathLossesWithGpsInfo{idxTrack}(:,[3,2]);
    curTxAzInDeg = TX_INFO_LOGS{1}(idxTrack).txAz;
    
    
    hpbwPolyshape = genBeamPolyshape(txLon, txLat, ...
        curTxAzInDeg, azimuthHpbwInDeg, BEAM_AREA_SIDE_LENGTH_IN_M);
    fnbwPolyshape = genBeamPolyshape(txLon, txLat, ...
        curTxAzInDeg, azimuthFnbwInDeg, BEAM_AREA_SIDE_LENGTH_IN_M);
    
    curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    colormap hot;
    plot3k([curRxLonLats, curPLs], ...
        'Marker', {'.', 12}, ...
        'ColorRange', EXPECTED_PATH_LOSS_RANGE, ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Track #', num2str(idxTrack)], ...
        'Longitude', 'Latitude', ...
        '', 'Path Loss (dB)'});
    xticks([]); yticks([]);
    
    plot(fnbwPolyshape, 'FaceColor', 'cyan');
    plot(hpbwPolyshape, 'FaceColor', 'green');
    
    axis(figAxisToSet); view(2);
    plotGoogleMapAfterPlot3k(curFig, 'satellite');
    
    curDirToSaveFig = fullfile(pathToSaveResults, ...
        ['basicPLsWithTxMainLobe_Track_', num2str(idxTrack), '.png']);
    saveas(curFig, curDirToSaveFig);
end

%% Time stamps on map.
startTimeInSForAllTracks = [1522520410 1522521163 1522522013 ...
    1522522690 1522523023 1522523676 ...
    1522524143 1522524548 1522524818 1522525036];
sampTimeStampsForAllTracks = cell(numOfTracks,1);
for idxT = 1:numOfTracks
    [curNumOfSamples, ~] = size(contiPathLossesWithGpsInfo{idxT});
    sampTimeStampsForAllTracks{idxT} = (startTimeInSForAllTracks(idxT) ...
        :(startTimeInSForAllTracks(idxT)+curNumOfSamples-1))';
end
allSampTimeStamps = vertcat(sampTimeStampsForAllTracks{:});
allRelSampTimeStamps = allSampTimeStamps-allSampTimeStamps(1);

curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
hTx = plot(TX_LON, TX_LAT, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);

colormap hot;
plot3k([allRxLonLats, allRelSampTimeStamps], ...
    'Marker', {'.', 12}, ...
    'CBLabels', numOfTicklabels, 'Labels', ...
    {'', 'Longitude', 'Latitude', ...
    '', 'Time Stamp (s)'});
xticks([]); yticks([]);

axis(figAxisToSet); view(2);
plotGoogleMapAfterPlot3k(curFig, 'satellite');

curDirToSaveFig = fullfile(pathToSaveResults, 'timeStampsOnMap');
saveas(curFig, [curDirToSaveFig, '.png']);
saveas(curFig, [curDirToSaveFig, '.fig']);

%% Signal Difference on Map
% We will draw the value:
%       Tx signal strength (dB) - Rx signal strength (dB)
% For each RX location on the map for debugging.

contiSignalDiffs = arrayfun(@(idx) ...
    contiPathLossesWithGpsInfo{idx}(:,1) ...
    - contiPathLossesExtraInfo.contiTxGains{idx} ...
    - contiPathLossesExtraInfo.contiRxGains{idx}, ...
    1:numOfTracks, 'UniformOutput', false);
allSigDiffs = vertcat(contiSignalDiffs{:});
EXPECTED_SIGNAL_DIFF_RANGE = [max(allSigDiffs), min(allSigDiffs)]; %[195, 50]; %

curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
hTx = plot(TX_LON, TX_LAT, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
colormap hot;
plot3k([allRxLonLats, allSigDiffs], ...
    'Marker', {'.', 12}, ...
    'ColorRange', EXPECTED_SIGNAL_DIFF_RANGE, ...
    'CBLabels', numOfTicklabels, 'Labels', ...
    {'Overview for Signal Difference', ...
    'Longitude', 'Latitude', ...
    '', 'Signal Difference (dB)'});
xticks([]); yticks([]);
axis(figAxisToSet); view(2);
plotGoogleMapAfterPlot3k(curFig, 'satellite');

curDirToSaveFig = fullfile(pathToSaveResults, ...
    ['Overview_SigDiffs.png']);
saveas(curFig, curDirToSaveFig);

for idxTrack = 1:numOfTracks
    curSigDiffs = contiPathLossesWithGpsInfo{idxTrack}(:,1);
    curRxLonLats = contiPathLossesWithGpsInfo{idxTrack}(:,[3,2]);
    curTxAzInDeg = TX_INFO_LOGS{1}(idxTrack).txAz;
    
    hpbwPolyshape = genBeamPolyshape(txLon, txLat, ...
        curTxAzInDeg, azimuthHpbwInDeg, BEAM_AREA_SIDE_LENGTH_IN_M);
    fnbwPolyshape = genBeamPolyshape(txLon, txLat, ...
        curTxAzInDeg, azimuthFnbwInDeg, BEAM_AREA_SIDE_LENGTH_IN_M);
    
    curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    colormap hot;
    plot3k([curRxLonLats, curSigDiffs], ...
        'Marker', {'.', 12}, ...
        'ColorRange', EXPECTED_SIGNAL_DIFF_RANGE, ...
        'CBLabels', numOfTicklabels, 'Labels', ...
        {['Track #', num2str(idxTrack)], ...
        'Longitude', 'Latitude', ...
        '', 'Signal Difference (dB)'});
    xticks([]); yticks([]);
    
    plot(fnbwPolyshape, 'FaceColor', 'cyan');
    plot(hpbwPolyshape, 'FaceColor', 'green');
    
    axis(figAxisToSet); view(2);
    plotGoogleMapAfterPlot3k(curFig, 'satellite');
    
    curDirToSaveFig = fullfile(pathToSaveResults, ...
        ['sigDiffsWithTxMainLobe_Track_', num2str(idxTrack), '.png']);
    saveas(curFig, curDirToSaveFig);
end

if FLAG_GEN_FIGS_SILENTLY
    close all
end

% EOF