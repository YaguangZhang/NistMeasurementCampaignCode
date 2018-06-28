% MANUALLYLOCATETREES Manually mark trees on
%
%   When available, parameters including series number, location, TxAz,
%   TxEl, RxAz, TxEl and note, will be loaded. We will also hardcode some
%   constant parameters here. For exmample, please make sure the path to
%   LiDAR data is correct specified by ABS_PATH_TO_NIST_LIDAR_LAS.
%
%   We will use all capitalized variable names for the parameters we
%   generate here.
%
% Yaguang Zhang, Purdue, 04/09/2018

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% The UTM zone from the .xml file.
UTM_ZONE = '13 N';

% The absolute path to the .las file.
switch getenv('computername')
    case 'ZYGLABS-DELL'
        % ZYG's Dell laptop.
        ABS_PATH_TO_NIST_LIDAR_LAS ...
            = 'C:\Users\Zyglabs\OneDrive - purdue.edu\EARS - Simulations\3D Models\USGS\NIST\CO_SoPlatteRiver_Lot5_2013_001049.las';
    case ''
        % Expected to be Lemma the Mac machine in ZYG's lab.
        assert(ismac, unknownComputerErrorMsg);
        ABS_PATH_TO_NIST_LIDAR_LAS = fullfile(absPathMacLemma, ...
            'CO_SoPlatteRiver_Lot5_2013_001048.las');
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% The absolute path to save elevation information for the lidar data.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'ManuallyLocateTrees');
% The tree location GPS records exported from an Android app.
ABS_PATH_TO_TREE_LOCS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'Data', 'NIST foliage analysis tree locations.csv');

% Set this to be true to also show the tree locations recorded by the
% Android app.
FLAG_SHOW_RECOREDED_TREE_LOCS = true;
pathToSaveTreeLocsRecorded = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'FoliageAttenuationEstimation', ...
    'treeLocs.mat');
% Set this to be true to try to fetch elevation data from Google Maps.
% Currently, we do not have enough quote for this.
FLAG_FTECH_ELE_FROM_GOOGLE = false;

% Configure other paths accordingly.
[~, NIST_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_NIST_LIDAR_LAS);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [NIST_LIDAR_LAS_FILENAME, '_Elevations.mat']);
ABS_PATH_TO_SAVE_TREE_LOCS = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'treeLocs.mat');

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

%% Load the Lidar Data

disp(' ')
disp('    Loading LiDAR data ...')
lidarData = lasdata(ABS_PATH_TO_NIST_LIDAR_LAS);

lidarXYZ = [lidarData.x, lidarData.y, lidarData.z];
% The stored intensity is unit16 and we need double for plotting.
lidarIntensity = double(lidarData.get_intensity);

lidarNumSamps = length(lidarData.z);
disp('    Done!')

%% Load the (lat, lon, ele) Information

disp(' ')
disp('    Generating elevation information ...')

if exist(ABS_PATH_TO_SAVE_ELEVATIONS, 'file')==2
    disp('        Loading previous results ...')
    load(ABS_PATH_TO_SAVE_ELEVATIONS);
else
    disp('        Converting UTM to (lat, lon) ...')
    % UTM to (lat, lon).
    [lidarLats, lidarLons] = utm2deg(lidarData.x, lidarData.y, ...
        repmat(UTM_ZONE, lidarNumSamps, 1));
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, 'lidarLats', 'lidarLons');
end

if FLAG_FTECH_ELE_FROM_GOOGLE ...
        && (~exist('lidarAlts', 'var') || isnan(lidarAlts))
    disp('        Fetching elevation information from Google Maps ...')
    % Fetch the elevation data from Google Maps.
    GOOGLE_MAPS_API = 'AIzaSyDlkaE_QJxvRJpTutWG0N-LCvoT0e7FPHE';
    lidarAlts = nan;
    countTrials = 0;
    while isnan(lidarAlts)
        countTrials = countTrials+1;
        disp(['            Trial # ', num2str(countTrials), ' ...']);
        lidarAlts = getElevationsFromGoogle(lidarLats, lidarLons, ...
            GOOGLE_MAPS_API);
    end
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, 'lidarAlts', '-append');
end

disp('    Done!')

%% Generate a Few Plots for the LiDAR Data

disp(' ')
disp('    Generating plots for LiDAR Data ...')

hIntensity = figure;
plot3k([lidarData.x, lidarData.y, lidarIntensity]);
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('Intensity');
title('LiDAR Intenstiy in 3D');

pathToSaveFig = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'lidarIntensity.png');
saveas(hIntensity, pathToSaveFig);

hZ = figure;
plot3k(lidarXYZ);
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('LiDAR Z in 3D');

pathToSaveFig = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'lidarZ.png');
saveas(hZ, pathToSaveFig);

downSampFactor = 20;
hZDownSamped = figure;
plot3k(lidarXYZ(1:downSampFactor:end, :));
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(['LiDAR Z (Down sampled by ', num2str(downSampFactor) ') in 3D']);

pathToSaveFig = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    ['lidarZDownSamped_', num2str(downSampFactor), '.png']);
saveas(hZDownSamped, pathToSaveFig);

numBins = 100;
hZHist = figure;
hist(lidarData.z, numBins);
grid minor; xlabel('z (m)'); ylabel('Frequency (#)');
title('Histogram for LiDAR Z');

pathToSaveFig = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'lidarZHist.png');
saveas(hZHist, pathToSaveFig);

disp('    Done!')

%% Interactive Figure for Labeling Tree Locations


disp(' ')
disp('    Generating interactive plot for marking trees ...')

% For filtering out outliers and samples we do not care.

zRangeToShow = [1660, 1920]; % Hard coded.
latRangeToShow = [39.989188, 39.992223];
lonRangeToShow = [-105.278014, -105.273414];

lidarLonLatZToShow = [lidarLons, lidarLats, lidarData.z];
lidarLonLatZToShow((lidarLonLatZToShow(:,1)<lonRangeToShow(1) ...
    | lidarLonLatZToShow(:,1)>lonRangeToShow(2)), :) = [];
lidarLonLatZToShow((lidarLonLatZToShow(:,2)<latRangeToShow(1) ...
    | lidarLonLatZToShow(:,2)>latRangeToShow(2)), :) = [];
lidarLonLatZToShow((lidarLonLatZToShow(:,3)<zRangeToShow(1) ...
    | lidarLonLatZToShow(:,3)>zRangeToShow(2)), :) = [];
% To better plot the map, we nomalize z to [0, 1].
normalizedZ = lidarLonLatZToShow(:,3);
normalizedZ = normalizedZ-min(normalizedZ);
normalizedZ = normalizedZ/max(normalizedZ);
lidarLonLatZToShow(:,3) = normalizedZ;

if exist('hInterTreeMarker', 'var') && isvalid(hInterTreeMarker)
    close(hInterTreeMarker);
end

hInterTreeMarker = figure('CloseRequestFcn',@saveMarkLocs, ...
    'visible','off');
hold on;
plot3k(lidarLonLatZToShow, 'Marker',{'.',1});
grid on; xlabel('Longitude'); ylabel('Latitude'); zlabel('Nomalized Z');
title('Selected LiDAR Z in 3D'); xticks([]); yticks([]);
axis([lonRangeToShow, latRangeToShow, 0, 1]);
plot_google_map('MapType', 'satellite');

pathToSaveFig = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'lidarZSelected.png');
saveas(hInterTreeMarker, pathToSaveFig);

disp('        Fetching tree locations ...')
if exist(ABS_PATH_TO_SAVE_TREE_LOCS, 'file')==2
    load(ABS_PATH_TO_SAVE_TREE_LOCS);
else
    markLocs = [];
    save(ABS_PATH_TO_SAVE_TREE_LOCS, 'markLocs');
end

if FLAG_SHOW_RECOREDED_TREE_LOCS
    treeLocsFromTablet = loadGpsMarkersWithAlt(ABS_PATH_TO_TREE_LOCS, 'Marker*', ...
        pathToSaveTreeLocsRecorded);
    [numTreeLocations, ~] = size(treeLocsFromTablet);
    hTreeLocationRecords = plot3(treeLocsFromTablet(:,2), treeLocsFromTablet(:,1), ...
        ones(numTreeLocations,1), ...
        'bo', 'LineWidth', 2);
end

% Overlay trees on top of the figure.
[numTrees, ~] = size(markLocs);
if numTrees>0
    hTreeLocs = plot3(markLocs(:, 2), markLocs(:, 1), ones(numTrees,1), ...
        'rx', 'LineWidth', 2);
else
    hTreeLocs = plot3([],[],[]);
end
% Change to X-Y view.
view(2);

% Add the interactive function.
hInteractiveArea = fill3([lonRangeToShow(1), lonRangeToShow(1), ...
    lonRangeToShow(2), lonRangeToShow(2)], ...
    [latRangeToShow, latRangeToShow(end:-1:1)], ...
    ones(1,4), 'r', 'FaceColor','none');
if FLAG_SHOW_RECOREDED_TREE_LOCS
    legend([hTreeLocs, hTreeLocationRecords], ...
        'Manual Labels', 'Android GPS');
end
set(hInteractiveArea, 'EdgeColor', 'r', 'LineWidth', 2);
set(hInteractiveArea,'ButtonDownFcn', @(src,evnt) ...
    updateTreeMarkerState(src, evnt), ...
    'PickableParts','all', 'HitTest','on');

% Show the figure.
set(hInterTreeMarker, 'visible','on');

disp('    The interactive tool for manually marking the tree locations is ready!')
disp('    Please (left) click on the plot to add new tree locations ...')
disp(' ')
disp('    It''s OK to zoom in the figure and move around if necessary. ')
disp('    It''s also OK to manually modify the tree locations stored in the variable markLocs in the base workspace. ')
disp('    The variable markLocs will eventually be saved into a .mat file when the figure is closed. ')
disp(' ')
disp('    Done!')

% EOF