% PLOTALLGPSSAMPSONMAP Plot all GPS samples logged by the USRP on a map.
%
% Yaguang Zhang, Purdue, 03/31/2018

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_FOLIAGE_DATA = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'Data', '20180331_NistFoliage');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults');

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

%% Read In the Log Files

disp(' ---------------------- ')
disp('  plotAllGpsSampsOnMap ')
disp(' ---------------------- ')
disp(' ')
disp('    Load GPS samples...')
% Scan the series folder for GPS log files.
allGpsLogFiles = rdir(fullfile(ABS_PATH_TO_FOLIAGE_DATA, '**', '*.log'));
allGpsSamps = arrayfun(@(l) parseGpsLog(l.name), allGpsLogFiles);
allGpsSamps = allGpsSamps(arrayfun(@(s) ~isempty(s.gpsLocation), allGpsSamps));

% Parse the GPS log.
numGpsSamps = length(allGpsSamps);
[allLats, allLons] = deal(nan(numGpsSamps, 1));
for idxS = 1:numGpsSamps
    gpsLog = allGpsSamps(idxS);
    [lat, lon, ~, ~] = parseNmeaStr(gpsLog.gpsLocation);
    
    allLats(idxS) = lat;
    allLons(idxS) = lon;
end

% Load the TX location.
[markerLats, markerLons, markerNames] ...
    = loadGpsMarkers(absPathToGpsMarkerCsvFile);
idxTxMarker = find(strcmp(markerNames, 'Tx'));
latTx = markerLats(idxTxMarker);
lonTx = markerLons(idxTxMarker);

disp('    Done!')

%% Plot the Overview

disp(' ')
disp('    Plotting overview...')

hFigOverview = figure;
hold on;
hTx = plot(lonTx, latTx, 'g^');
hSamps = plot(allLons, allLats, 'r.');
axis tight;
plot_google_map('MapType', 'satellite');
legend([hTx, hSamps], 'Tx', 'Measurements');

% Save the plot.
pathToSaveOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Overview.png');
saveas(hFigOverview, pathToSaveOverviewPlot);

disp('    Done!')

disp('    Done!')

%% Plot Each Continuous Recording

disp(' ')
disp('    Plotting separate recording activities...')

pathToSaveSeparateRecs = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'GpsSampsOnMap');

% Find all GnuRadio .out files.
allSigOutFiles = rdir(fullfile(ABS_PATH_TO_FOLIAGE_DATA, '**', '*.out'));
% We do not need the filtered version.
allSigOutFiles = allSigOutFiles(arrayfun(@(p) ...
    ~isempty(regexp(p.name, '\d+.out','match')), allSigOutFiles));
numAllSigOutFiles = length(allSigOutFiles);
for idxRec = 1:numAllSigOutFiles
    curOutFilePath = allSigOutFiles(idxRec).name;    
    curGpsLogs = rdir(fullfile(allSigOutFiles(idxRec).folder, '*.log'));

    % Fetch the (lat, lon) from these GPS logs.
    [curLats, curLons, ~, ~] = parseGpsLogs(curGpsLogs);

    % Plot.
    hContiRec = figure;
    hold on;
    hTx = plot(lonTx, latTx, 'g^');
    hSamps = plot(allLons, allLats, 'r.');
    hCurRec = plot(curLons, curLats, 'b.');
    axis tight;
    plot_google_map('MapType', 'satellite');
    legend([hTx, hCurRec, hSamps], ...
        'Tx', ['Route #', num2str(idxRec)], 'Other measurements');

    % Save the plot.
    pathToSaveCurRecPlot = fullfile(pathToSaveSeparateRecs, ...
        ['Overview_', num2str(idxRec), '.png']);
    saveas(hContiRec,  pathToSaveCurRecPlot);
end

disp('    Done!')

% EOF