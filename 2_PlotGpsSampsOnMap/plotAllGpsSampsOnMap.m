% PLOTALLGPSSAMPSONMAP Plot all GPS samples logged by the USRP on a map.
%
% Yaguang Zhang, Purdue, 03/31/2018

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_FOLIAGE_DATA = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'Data', '20180331_NistFoliage');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
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
[allLats, allLons, allAlts] = deal(nan(numGpsSamps, 1));
for idxS = 1:numGpsSamps
    gpsLog = allGpsSamps(idxS);
    [lat, lon, alt, ~] = parseNmeaStr(gpsLog.gpsLocation);
    
    allLats(idxS) = lat;
    allLons(idxS) = lon;
    allAlts(idxS) = alt;
end

% Load the TX location.
[markerLats, markerLons, markerNames] ...
    = loadGpsMarkers(absPathToGpsMarkerCsvFile);
idxTxMarker = find(strcmp(markerNames, 'Tx'));
latTx = markerLats(idxTxMarker);
lonTx = markerLons(idxTxMarker);

disp('    Done!')

%% Plot 2D Overview

disp(' ')
disp('    Plotting 2D overview...')

hFigOverview = figure;
hold on;
hTx = plot(lonTx, latTx, 'g^');
hSamps = plot(allLons, allLats, 'r.');
axis tight;
plot_google_map('MapType', 'satellite');
legend([hTx, hSamps], 'Tx', 'Measurements');
title('GPS Sample Overview on Map');

% Save the plot.
pathToSaveOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Overview.png');
saveas(hFigOverview, pathToSaveOverviewPlot);

disp('    Done!')

%% Plot 3D Overview

disp(' ')
disp('    Plotting 3D overview...')

altTxGoogle = 1778.9871826;
altRxMin = min(allAlts);

hFigOverview3D = figure;
hold on;
% Two possible TX locaitons with different altitudes.
hTxGoogle = plot3(lonTx, latTx, altTxGoogle, 'g^');
hRxMin = plot3(lonTx, latTx, altRxMin, 'b^');
hSamps = plot3(allLons, allLats, allAlts, 'r.');
legend([hTxGoogle, hRxMin, hSamps], ...
    'Tx By Google Alt', 'Tx by Rx Min Alt', 'Measurements');
title('GPS Sample Overview in 3D'); view(15, 30); grid minor;
xlabel('x (m)'); ylabel('y (m)'); zlabel('Altitude (m)');

% Save the plot.
pathToSaveOverview3DPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Overview3D.png');
saveas(hFigOverview3D, pathToSaveOverview3DPlot);

disp('    Done!')

%% Plot Each Continuous Recording

disp(' ')
disp('    Plotting separate recording activities...')

pathToSaveSeparateRecs = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'GpsSampsOnMap');
% Create directories if necessary.
if exist(pathToSaveSeparateRecs, 'dir')~=7
    mkdir(pathToSaveSeparateRecs);
end

% Find all GnuRadio .out files.
allSigOutFiles = rdir(fullfile(ABS_PATH_TO_FOLIAGE_DATA, '**', '*.out'));
% We do not need the filtered version.
allSigOutFiles = allSigOutFiles(arrayfun(@(p) ...
    ~isempty(regexp(p.name, '\d+.out','match')), allSigOutFiles));
numAllSigOutFiles = length(allSigOutFiles);
for idxRec = 1:numAllSigOutFiles
    curSeriesFolderPath = allSigOutFiles(idxRec).folder;
    
    % Get the Serires numbet just to make sure.
    matchedSeriesNums = regexp(curSeriesFolderPath, 'Series_(\d+)', 'tokens');
    seriesNum = str2double(matchedSeriesNums{1}{1});
    
    % Find all the GPS logs in the current folder.
    curOutFilePath = allSigOutFiles(idxRec).name;    
    curGpsLogs = rdir(fullfile(curSeriesFolderPath, '*.log'));

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
        'Tx', ['Route #', num2str(seriesNum)], 'Other measurements');

    % Save the plot.
    pathToSaveCurRecPlot = fullfile(pathToSaveSeparateRecs, ...
        ['Overview_', num2str(seriesNum), '.png']);
    saveas(hContiRec,  pathToSaveCurRecPlot);
end

disp('    Done!')

% EOF