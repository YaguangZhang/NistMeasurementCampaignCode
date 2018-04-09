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
ABS_PATH_TO_DATA = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
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
allGpsLogFiles = rdir(fullfile(ABS_PATH_TO_DATA, '**', '*.log'));
allGpsSamps = arrayfun(@(l) parseGpsLog(l.name), allGpsLogFiles);
allGpsSamps = allGpsSamps(arrayfun(@(s) ~isempty(s.gpsLocation), allGpsSamps));

% Parse the GPS log.
numGpsSamps = length(allGpsSamps);
[lats, lons] = deal(nan(numGpsSamps, 1));
for idxS = 1:numGpsSamps
    gpsLog = allGpsSamps(idxS);
    [lat, lon, ~, ~] = parseNmeaStr(gpsLog.gpsLocation);
    
    lats(idxS) = lat;
    lons(idxS) = lon;
end
disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

hFigOverview = figure;
plot(lons, lats, 'r.');
plot_google_map('MapType', 'satellite');

% Save the plot.
pathToSaveOverviewPlot = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Overview.png');
saveas(hFigOverview, pathToSaveOverviewPlot);

disp('    Done!')
% EOF