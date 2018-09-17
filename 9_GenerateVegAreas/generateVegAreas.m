% GENERATEVEGAREAS Generate vegetation areas on a map for the NIST data
% set.
%
%   When available, parameters including series number, location, TxAz,
%   TxEl, RxAz, TxEl and note, will be loaded. We will also hardcode some
%   constant parameters here. For exmample, please make sure the path to
%   LiDAR data is correct specified by ABS_PATH_TO_NIST_LIDAR_LAS.
%
%   We will use all capitalized variable names for these parameters.
%
%   The resulting data for the vegetation areas, i.e. vegAreas, will
%   essentially be for an image, which can be view by imtool(vegAreas),
%   with vegetation areas highlighed in a black background.
%
% Yaguang Zhang, Purdue, 08/06/2018

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% The absolute path to the .las file.
switch getenv('computername')
    case 'ZYGLABS-DELL'
        % ZYG's Dell laptop.
        ABS_PATH_TO_NIST_LIDAR_LAS ...
            = 'C:\Users\Zyglabs\OneDrive - purdue.edu\EARS - Simulations\3D Models\USGS\NIST\CO_SoPlatteRiver_Lot5_2013_001049.las';
    case 'ARTSY'
        % ZYG's lab desktop.
        ABS_PATH_TO_NIST_LIDAR_LAS ...
            = 'C:\Users\YaguangZhang\OneDrive - purdue.edu\EARS - Simulations\3D Models\USGS\NIST\CO_SoPlatteRiver_Lot5_2013_001049.las';
    case ''
        % Expected to be Lemma the Mac machine in ZYG's lab.
        assert(ismac, unknownComputerErrorMsg);
        ABS_PATH_TO_NIST_LIDAR_LAS = fullfile(absPathMacLemma, ...
            'CO_SoPlatteRiver_Lot5_2013_001048.las');
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% The absolute path to save results.
ABS_PATH_TO_TREE_MARKER_RES = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'ManuallyLocateTrees');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'AutoGeneratedVegAreas');

% Configure other paths accordingly.
[~, NIST_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_NIST_LIDAR_LAS);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(ABS_PATH_TO_TREE_MARKER_RES, ...
    [NIST_LIDAR_LAS_FILENAME, '_Elevations.mat']);
ABS_PATH_TO_SAVE_VEG_AREAS_IMG = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'vegAreas.png');
ABS_PATH_TO_SAVE_VEG_AREAS_META = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'vegAreasMeta.mat');

% The UTM zone from the .xml file is 13N, which indicates the zone 13 on
% the Northern Hemisphere. However, we need the latitudinal band (S), too.
UTM_ZONE = '13 S'; % The UTM zone expected for the NIST data.

% For filtering out outliers and samples we do not care.
Z_RANGE = [1660, 1920]; % Hard coded.
LAT_RANGE = [39.989188, 39.992223];
LON_RANGE = [-105.278014, -105.273414];

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

%% Get GPS Locations of the Lidar Data and Set Parameters

disp(' ')
disp('    Loading results from: ')
disp('      - xxx_Elevations.mat')

assert(exist(ABS_PATH_TO_SAVE_ELEVATIONS, 'file')==2, ...
    'Couldn''t find xxx_Elevations.mat! Please run PostProcessing/7_ManuallyLocateTrees/manuallyLocateTrees.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get lidarLats, lidarLons, nistElevData, and if available, lidarAlts.
load(ABS_PATH_TO_SAVE_ELEVATIONS);

%% Interactive Figure for Labeling Tree Locations

disp(' ')
disp('    Generating interactive plot for drawing vegetation areas ...')

% We will set the image resolution according to LAT_RANGE and LON_RANGE to
% get roughly 0.1x0.1-square-meter pixels.
pixelRecWidth = 0.1; % In meters.
VEG_AREA_IMG_RESOLUTION = [nan nan]; % (Width, Height) in pixels.
% Convert (lon, lat) boundary points to UTM.
[xBottomLeftBound, yBottomLeftBound, zoneBottomLeftBound] ...
    = deg2utm(min(LAT_RANGE), min(LON_RANGE));
[xTopRightBound, yTopRightBound, zoneTopRightBound] ...
    = deg2utm(max(LAT_RANGE), max(LON_RANGE));
assert(strcmp(zoneBottomLeftBound, UTM_ZONE) ...
    && strcmp(zoneTopRightBound, UTM_ZONE), ...
    ['Boundary is expected to be in the UTM zone ', UTM_ZONE, '!']);

% Image width.
VEG_AREA_IMG_RESOLUTION(1) = ceil(...
    (xTopRightBound-xBottomLeftBound)./pixelRecWidth +1);
% Image height.
VEG_AREA_IMG_RESOLUTION(2) = ceil( ...
    (yTopRightBound-yBottomLeftBound)./pixelRecWidth +1);

% Limit the lidar data to the range we care about.
lidarLonLatAltZToShow = [lidarLons, lidarLats, lidarAlts, lidarData.z];
lidarLonLatAltZToShow((lidarLonLatAltZToShow(:,1)<LON_RANGE(1) ...
    | lidarLonLatAltZToShow(:,1)>LON_RANGE(2)), :) = [];
lidarLonLatAltZToShow((lidarLonLatAltZToShow(:,2)<LAT_RANGE(1) ...
    | lidarLonLatAltZToShow(:,2)>LAT_RANGE(2)), :) = [];
lidarLonLatAltZToShow((lidarLonLatAltZToShow(:,3)<Z_RANGE(1) ...
    | lidarLonLatAltZToShow(:,3)>Z_RANGE(2)), :) = [];

% Plot z on a map. To better show the data, we nomalize z to [0, 1].
normalizedZ = lidarLonLatAltZToShow(:,4);
normalizedZ = normalizedZ-min(normalizedZ);
normalizedZ = normalizedZ/max(normalizedZ);
lidarLonLatNormZToShow = [lidarLonLatAltZToShow(:,1:2) normalizedZ];

if exist('hInterVegAreaMarker', 'var') && isvalid(hVegAreaOnMap)
    close(hVegAreaOnMap);
end

hVegAreaOnMap = figure('visible','off');
hold on;
plot3k(lidarLonLatNormZToShow, 'Marker',{'.',1});
grid on; xlabel('Longitude'); ylabel('Latitude'); zlabel('Nomalized Z');
title('Vegetation Areas with Colored LiDAR z Values (Normalized)'); 
xticks([]); yticks([]);
axis([LON_RANGE, LAT_RANGE, 0, 1]);
plot_google_map('MapType', 'satellite');

% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hVegAreaOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

disp('        Fetching vegetation areas ...')
if exist(ABS_PATH_TO_SAVE_VEG_AREAS_META, 'file')==2
    disp('             Loading history restuls ...')
    load(ABS_PATH_TO_SAVE_VEG_AREAS_META);
else
    disp('             Initializing meta data ...')
    
    % We will generate memoization tables for quickly fetch the GPS
    % locations and UTM coordinates for the pixels. The information is
    % stored in a structure VEG_AREA_IMG_META with fields:
    %   - LATS, LONS
    %     The GPS positions arranged in 2D image formats.
    %   - XS, YS
    %     The UTM coordinates arranged in 2D image formats.
    %   - UTM_ZONE
    %     A string specifying the UTM zone.
    %   - ALTS
    %     The bare earth altitudes arranged in a 2D image format, which are
    %     estimated by interpolating the Goolge elevation results for the
    %     Lidar data.
    %   - ZS
    %     The Lidar Z values arranged in a 2D image format, which are
    %     estimated by interpolating the NIST Lidar data.
    VEG_AREA_IMG_META.IMG_RESOLUTION = VEG_AREA_IMG_RESOLUTION;
    VEG_AREA_IMG_META.LAT_RANGE = LAT_RANGE;
    VEG_AREA_IMG_META.LON_RANGE = LON_RANGE;
    
    [VEG_AREA_IMG_META.LATS, VEG_AREA_IMG_META.LONS, ...
        VEG_AREA_IMG_META.XS, VEG_AREA_IMG_META.YS] ...
        = deal(nan(VEG_AREA_IMG_RESOLUTION([2,1])));
    VEG_AREA_IMG_META.UTM_ZONE = UTM_ZONE;
    % Build the grid in UTM space and compute the cooresponding GPS lats
    % and lons.
    [vegAreaImgMetaLatsCache, vegAreaImgMetaLonsCache, ...
        vegAreaImgMetaXsCache, vegAreaImgMetaYsCache ] ...
        = deal(nan(VEG_AREA_IMG_RESOLUTION([2,1])));
    disp('                 Computing position coordinates for the grid ...')
    parfor idxX = 1:VEG_AREA_IMG_RESOLUTION(1)
        % Progress feedback for every percent of work.
        normProgPerc = idxX./VEG_AREA_IMG_RESOLUTION(1).*100;  %#ok<PFBNS>
        if rem(idxX, VEG_AREA_IMG_RESOLUTION(1)/100) == 0
            disp(['                     ', num2str(normProgPerc), '% ...'])
        end
        [curImgLats, curImgLons] ...
            = pixIndices2LatLon( ...
            ones(VEG_AREA_IMG_RESOLUTION(2),1).*idxX, ...
            (1:VEG_AREA_IMG_RESOLUTION(2))', ...
            LAT_RANGE, LON_RANGE, VEG_AREA_IMG_RESOLUTION);
        
        [vegAreaImgMetaXsCache(:, idxX), ...
            vegAreaImgMetaYsCache(:, idxX), curUtmZones] ...
            = deg2utm(curImgLats, curImgLons);
        
        assert(all(arrayfun(@(idx) ...
            strcmp(curUtmZones(idx,:), UTM_ZONE), ...
            1:length(curImgLats))), ...
            'All the locations should be in the same UTM zone!');
        
        vegAreaImgMetaLatsCache(:, idxX) ...
            = curImgLats;
        vegAreaImgMetaLonsCache(:, idxX) ...
            = curImgLons;
    end
    VEG_AREA_IMG_META.XS = vegAreaImgMetaXsCache;
    VEG_AREA_IMG_META.YS = vegAreaImgMetaYsCache;
    VEG_AREA_IMG_META.LATS = vegAreaImgMetaLatsCache;
    VEG_AREA_IMG_META.LONS = vegAreaImgMetaLonsCache;
    
    disp('                 Interperlating data for lats ...')
    % Interpolate the Lidar data lidarLonLatAltZToShow to a grid
    % corresponding to the image.
    [gridPixXs, gridPixYs] = meshgrid(1:VEG_AREA_IMG_RESOLUTION(1), ...
        1:VEG_AREA_IMG_RESOLUTION(2));
    % Interperlating for alts.
    [nistElevDataLons, nistElevDataLats] = meshgrid( ...
        nistElevData.longs, nistElevData.lats);
    % Note: Grid lats, lons, xs and ys are essentially the same as those
    % for VEG_AREA_IMG_META.
    VEG_AREA_IMG_META.ALTS = interp2( ... % Lons, lats, alts.
        nistElevDataLons, nistElevDataLats, ...
        nistElevData.elev, ...
        VEG_AREA_IMG_META.LONS, VEG_AREA_IMG_META.LATS);
    
    disp('                 Interperlating data for zs ...')
    % Interperlating for zs.
    getInterZ = scatteredInterpolant( ...
        lidarLonLatAltZToShow(:,1), ...
        lidarLonLatAltZToShow(:,2), ...
        lidarLonLatAltZToShow(:,4));
    VEG_AREA_IMG_META.ZS = getInterZ( ...
        VEG_AREA_IMG_META.LONS, VEG_AREA_IMG_META.LATS);
    
    % Initialize vegAreas.
    treeHeight = VEG_AREA_IMG_META.ZS - VEG_AREA_IMG_META.ALTS;
    
    disp('                 Generating figures for different minTreeHs ...')
    % We will generate .png figures for different tree height threshold to
    % choose the best one.
    minTreeHs = [0:0.1:3 3.5 4:7 9 11]; % In meter.
    minTreeHFigHeight = 1500; % In Pixel.
    for curMinTreeH = minTreeHs
        hMinTreeHFig = figure('Position', ...
            [0 minTreeHFigHeight ...
            ceil(1500./VEG_AREA_IMG_RESOLUTION(2) ...
            .*VEG_AREA_IMG_RESOLUTION(1)) ...
            minTreeHFigHeight], ...
            'visible','off'); hold on;
        boolsIsVegArea = treeHeight(:)>=curMinTreeH;
        plot3k([VEG_AREA_IMG_META.LONS(boolsIsVegArea), ...
            VEG_AREA_IMG_META.LATS(boolsIsVegArea), ...
            treeHeight(boolsIsVegArea)])
        axis([LON_RANGE, LAT_RANGE, 0 max(treeHeight(:))]);
        FLAG_GOOGLE_MAP_FETCHED = false;
        while ~FLAG_GOOGLE_MAP_FETCHED
            try
            plot_google_map('Maptype', 'satellite'); 
            FLAG_GOOGLE_MAP_FETCHED = true;
            catch
                disp('Error in overlaying Google map! Retrying ... ');
            end
        end
        view(2);
        % The command plot_google_map messes up the color legend of plot3k,
        % so we will have to fix it here.
        hCb = findall( allchild(hMinTreeHFig), 'type', 'colorbar');
        hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
        xlabel('Longitude'); ylabel('Latitude'); xticks([]); yticks([]);
        export_fig(fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['vegArea_minTreeHeight_', ...
            num2str(curMinTreeH, '%.1f'), '_m.png']), ...
            '-png', '-transparent', hMinTreeHFig);
        close(hMinTreeHFig);
    end
    
    % We will save the vegAreas as a plain 2D image (width x height), with
    % the background being black (0) and vegetable area being white (1).
    vegAreas = zeros(VEG_AREA_IMG_RESOLUTION([2,1]));
    % The RX height is ~0.5m.
    VEG_AREA_IMG_META.MIN_TREE_H = 0.5;
    boolsIsVegArea = treeHeight(:)>=VEG_AREA_IMG_META.MIN_TREE_H;
    vegAreas(boolsIsVegArea) = 1;
    
    save(ABS_PATH_TO_SAVE_VEG_AREAS_META, 'vegAreas', 'LAT_RANGE', ...
        'LON_RANGE', 'VEG_AREA_IMG_RESOLUTION', 'VEG_AREA_IMG_META', ...
        '-v7.3');
end

disp('        Plotting ...')
% Color the vegetation areas on the figure.
hVegAreas = plotVegAreasOnMap(hVegAreaOnMap, vegAreas, ...
    VEG_AREA_IMG_META, 1);

% Change to X-Y view.
view(2);

% Add a square boundary.
hAreaOfInterest = fill3([LON_RANGE(1), LON_RANGE(1), ...
    LON_RANGE(2), LON_RANGE(2)], ...
    [LAT_RANGE, LAT_RANGE(end:-1:1)], ...
    ones(1,4), 'r', 'FaceColor','none');
set(hAreaOfInterest, 'EdgeColor', ones(3,1).*0.9, 'LineWidth', 2);

% Show the figure.
set(hVegAreaOnMap, 'visible','on');
drawnow;
legend(hVegAreas, 'Vegetation area');
transparentizeCurLegends;

% Save the illustration figure.
saveas(hVegAreaOnMap, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, 'vegAreaOverview.png'));

disp('    Done!')