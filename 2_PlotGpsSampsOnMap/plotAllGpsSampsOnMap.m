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
    'PostProcessingResults', 'AllGpsSampsOnMap');

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

%% Read In the Log Files

% Scan the series folder for GPS log files.
allGpsLogFiles = rdir(fullfile(ABS_PATH_TO_DATA, '**', '*.log'));
allGpsSamps = arrayfun(@(l) parseGpsLog(l.name), allGpsLogFiles);
allGpsSamps = allGpsSamps(arrayfun(@(s) ~isempty(s.gpsLocation), allGpsSamps));

% Parse the GPS log.
numGpsSamps = length(allGpsSamps);
[lats, lons] = deal(nan(numGpsSamps, 1));
for idxS = 1:numGpsSamps
    gpsLog = allGpsSamps(idxS);
    gpsLogSample = nmealineread(gpsLog.gpsLocation);
    
    lat = gpsLogSample.latitude;
    lon = gpsLogSample.longitude;
    alt = gpsLogSample.altitude;
    
    % Add a minus sign if it is W or S.
    if(isW(gpsLog.gpsLocation))
        lon = -lon;
    end
    if(isS(gpsLog.gpsLocation))
        lat = -lat;
    end
    
    lats(idxS) = lat;
    lons(idxS) = lon;
end

%% Plot
plot(lons, lats, 'r.');
plot_google_map('MapType', 'satellite');

% EOF