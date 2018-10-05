% MANUALLYLOCATETREES Manually mark trees on a map for the NIST data set.
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

% The UTM zone from the .xml file is 13N, which indicates the zone 13 on
% the Northern Hemisphere. However, we need the latitudinal band (S), too.
UTM_ZONE = '13 S';

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

% The absolute path to save elevation information for the lidar data.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'ManuallyLocateTrees');
% The tree location GPS records exported from an Android app.
ABS_PATH_TO_TREE_LOCS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'Data', 'NIST foliage analysis tree locations.csv');

% Set this to be true to also show the tree locations recorded by the
% Android app.
FLAG_SHOW_RECOREDED_TREE_LOCS = false;
pathToSaveTreeLocsRecorded = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'FoliageAttenuationEstimation', ...
    'treeLocs.mat');
% Set this to be true to try to fetch elevation data from Google Maps.
% Currently, we do not have enough quota for this.
FLAG_FTECH_ELE_FROM_GOOGLE = false;

% Configure other paths accordingly.
[~, NIST_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_NIST_LIDAR_LAS);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [NIST_LIDAR_LAS_FILENAME, '_Elevations.mat']);
ABS_PATH_TO_SAVE_TREE_LOCS = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'treeLocs.mat');

% Set this flag to be true to disable the interactive function of the
% collector and generate overview plots instead.
FLAG_GENERATE_OVERVIEW_PLOTS_ONLY = true;

% For filtering out outliers and samples of the Lidar data we do not care.
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
    
    % Use USGS 1/3 arc-second (~10m) resolution data for US terrain
    % elevation.
    region = fetchregion(LAT_RANGE, LON_RANGE, 'display', true);
    nistElevData = region.readelevation(LAT_RANGE, LON_RANGE, ...
        'sampleFactor', 1, ...
        'display', true);
    % One can preview the elevation data fetched using:
    %       dispelev(nistElevData, 'mode', 'latlong');
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, 'lidarLats', 'lidarLons', ...
        'nistElevData');
end

% we will only try to generate the elevation information when it is not
% available.
if (~exist('lidarAlts', 'var') || all(isnan(lidarAlts)))
    if FLAG_FTECH_ELE_FROM_GOOGLE
        lidarAlts = nan;
        % Use Google as the elevation source for better resolution.
        % However, this approach is limited by the quota Google provides.
        disp('        Fetching elevation information from Google Maps ...')
        % Fetch the elevation data from Google Maps.
        GOOGLE_MAPS_API = 'AIzaSyDlkaE_QJxvRJpTutWG0N-LCvoT0e7FPHE';
        
        countTrials = 0;
        while isnan(lidarAlts)
            countTrials = countTrials+1;
            disp(['            Trial # ', num2str(countTrials), ' ...']);
            lidarAlts = getElevationsFromGoogle(lidarLats, lidarLons, ...
                GOOGLE_MAPS_API);
        end
    else
        % Interperlate nistElevData to get lidarAlts.
        [nistElevDataLons, nistElevDataLats] = meshgrid( ...
            nistElevData.longs, nistElevData.lats);
        lidarAlts = interp2( ...
            nistElevDataLons, nistElevDataLats, ...
            nistElevData.elev, ...
            lidarLons, lidarLats);
    end
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, 'FLAG_FTECH_ELE_FROM_GOOGLE', ...
        'lidarAlts', '-append');
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

if exist('hInterTreeMarker', 'var') && isvalid(hInterTreeMarker)
    close(hInterTreeMarker);
end

if FLAG_GENERATE_OVERVIEW_PLOTS_ONLY
    hInterTreeMarker = figure('visible','off');
else
    hInterTreeMarker = figure('CloseRequestFcn',@saveMarkLocs, ...
        'visible','off');
end
hold on;
plot3k(lidarLonLatNormZToShow, 'Marker',{'.',1});
grid on; xlabel('Longitude'); ylabel('Latitude'); zlabel('Nomalized Z');
title('Tree Locations with Colored LiDAR z Values (Normalized)');
xticks([]); yticks([]);
axis([LON_RANGE, LAT_RANGE, 0, 1]);
plot_google_map('MapType', 'satellite');

% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hInterTreeMarker), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

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
hInteractiveArea = fill3([LON_RANGE(1), LON_RANGE(1), ...
    LON_RANGE(2), LON_RANGE(2)], ...
    [LAT_RANGE, LAT_RANGE(end:-1:1)], ...
    ones(1,4), 'r', 'FaceColor','none');
if FLAG_SHOW_RECOREDED_TREE_LOCS
    legend([hTreeLocs, hTreeLocationRecords], ...
        'Manual Labels', 'Android GPS', 'Location', 'southeast', ...
        'AutoUpdate', 'off');
else
    legend(hTreeLocs, 'Manual Labels', 'Location', 'southeast', ...
        'AutoUpdate', 'off');
end
set(hInteractiveArea, 'EdgeColor', ones(3,1).*0.9, 'LineWidth', 2);
uistack(hTreeLocs,'top');
if ~FLAG_GENERATE_OVERVIEW_PLOTS_ONLY
    set(hInteractiveArea,'ButtonDownFcn', @(src,evnt) ...
        updateTreeMarkerState(src, evnt), ...
        'PickableParts','all', 'HitTest','on');
    
    disp('    The interactive tool for manually marking the tree locations is ready!')
    disp('    Please (left) click on the plot to add new tree locations ...')
    disp(' ')
    disp('    It''s OK to zoom in the figure and move around if necessary. ')
    disp('    It''s also OK to manually modify the tree locations stored in the variable markLocs in the base workspace. ')
    disp('    The variable markLocs will eventually be saved into a .mat file when the figure is closed. ')
    disp(' ')
    disp('    Done!')
end

% Show the figure.
set(hInterTreeMarker, 'visible','on');

if FLAG_GENERATE_OVERVIEW_PLOTS_ONLY
    disp('    Generating overview plots ... ')
    disp('        (Please set FLAG_GENERATE_OVERVIEW_PLOTS_ONLY to be ')
    disp('        false if you would like to use the interactive tree ')
    disp('        marker tool instead)')
    
    axisValuesForOverviews = { ...
        [-105.27824323, -105.27319200, 39.98900606, 39.99241554]; ...
        [-105.27701930, -105.27483020, 39.99008676, 39.99156435]; ...
        [-105.27661236, -105.27555726, 39.99085158, 39.99156374]};
    
    % For plotting the areas indicated by the next axis vector.
    hNextArea = nan;
    numOverviewAreas = length(axisValuesForOverviews);
    for idxOverview = 1:numOverviewAreas
        if isgraphics(hNextArea) && isvalid(hNextArea)
            delete(hNextArea);
        end
        
        curFileName = ['Overview_', num2str(idxOverview)];
        pathToSaveOverviews = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            curFileName);
        
        axis(axisValuesForOverviews{idxOverview});
        saveas(hInterTreeMarker, [pathToSaveOverviews, '.png']);
        
        % Export an .eps copy for papers.
        pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
            'PostProcessingResults', '1_EpsFigs');
        saveEpsFigForPaper(hInterTreeMarker, ...
            fullfile(pathToSavePaperFigs, ['2_1_', curFileName, '.eps']));
        
        if idxOverview<numOverviewAreas
            hNextArea = plot3(...
                axisValuesForOverviews{idxOverview+1}([1 1 2 2 1]), ...
                axisValuesForOverviews{idxOverview+1}([3 4 4 3 3]), ...
                ones(1,5), '-.', 'LineWidth', 1, 'Color', ones(3,1));
            saveas(hInterTreeMarker, [pathToSaveOverviews, '_WithNextArea.png']);
            
            % Export an .eps copy for papers.
            curFileName = [curFileName, '_WithNextArea'];
            saveEpsFigForPaper(hInterTreeMarker, ...
                fullfile(pathToSavePaperFigs, ['2_1_', curFileName, '.eps']));
        end
    end
    
    % Generate a large image with only satellite image and lidar data for
    % manually coloring leave area, which is used to calculate
    % site-specific vagetation depth for applying the ITU model.
    VEG_AREA_IMG_RESOLUTION = [2000 2000]; % In pixel.
    VEG_AREA_IMG_AXIS = [-105.27824323, -105.27319200, ...
        39.98900606, 39.99241554];
    
    figure(hInterTreeMarker);
    curInterTreeMarkerFigPos = get(hInterTreeMarker, 'Position');
    curInterTreeMarkerAxis = axis;
    set(hInterTreeMarker, 'Position', [0 0 VEG_AREA_IMG_RESOLUTION]);
    axis(VEG_AREA_IMG_AXIS);
    plot_google_map('MapType', 'satellite');
    vegAreaFrame = getframe(gcf);
    
    vegAreaImgPathToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        'vegArea.png');
    vegAreaImgParaPathToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        'vegAreaParameters.mat');
    
    imwrite(vegAreaFrame.cdata, vegAreaImgPathToSave);
    save(vegAreaImgParaPathToSave, 'VEG_AREA_IMG_RESOLUTION', ...
        'VEG_AREA_IMG_AXIS', 'vegAreaFrame');
    
    figure(hInterTreeMarker);
    set(hInterTreeMarker, 'Position', curInterTreeMarkerFigPos);
    axis(curInterTreeMarkerAxis);
    
    disp('    Overview plots have been generated ')
    disp(' ')
end

%% Generate an Tree Location and Foliage Area Overview in One Plot

% This is for the IEEE paper.
ABS_PATH_TO_VEG_AREAS_META = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'AutoGeneratedVegAreas', 'vegAreasMeta.mat');
if exist(ABS_PATH_TO_VEG_AREAS_META, 'file')==2
    disp(' ')
    disp('    Generating overview plot for all site-specific features ...')

    % We need the foliage areas marked by:
    %   9_GenerateVegAreas/generateVegAreas.m
    load(ABS_PATH_TO_VEG_AREAS_META);
    
    axisToShow = [-105.27701930, -105.27483020, 39.99008676, 39.99156435];
    
    % The left side will show the foliage areas and the right part will
    % have the lidar data with manually labeled tree markers.
    centerLon = mean(axisToShow(1:2));
    centerLat = mean(axisToShow(3:4));
    
    hVegOverview = figure('visible','off'); hold on;
    lidarLonLatNormZToShowRightHalf ...
        = lidarLonLatNormZToShow(lidarLonLatNormZToShow(:,1)>centerLon, :);
    
    normalizedZRightHalf = lidarLonLatNormZToShowRightHalf(:,3);
    normalizedZRightHalf = normalizedZRightHalf-min(normalizedZRightHalf);
    normalizedZRightHalf = normalizedZRightHalf/max(normalizedZRightHalf);
    lidarLonLatNormZToShowRightHalfNorm ...
        = [lidarLonLatNormZToShowRightHalf(:,1:2) ...
        normalizedZRightHalf];
    plot3k(lidarLonLatNormZToShowRightHalfNorm, 'Marker',{'.',1});
    grid on; xlabel('Longitude'); ylabel('Latitude'); zlabel('Nomalized Z');
    title('Foliage Areas and Tree Locations (with Normalized LiDAR Values)');
    xticks([]); yticks([]);
    axis(axisToShow);
    plot_google_map('MapType', 'satellite');
    
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = findall( allchild(hVegOverview), 'type', 'colorbar');
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    
    % Overlay trees on top of the figure.
    [numTrees, ~] = size(markLocs);
    boolsToShowMarkLocs ...
        = (markLocs(:,2)>centerLon)&(markLocs(:,1)<centerLat);
    hTreeLocs = plot3(markLocs(boolsToShowMarkLocs, 2), ...
        markLocs(boolsToShowMarkLocs, 1), ...
        ones(sum(boolsToShowMarkLocs),1), ...
        'rx', 'LineWidth', 1.5);
    
    % Center vertical divider.
    deltaAxisLonToShow = axisToShow(2)-axisToShow(1);
    deltaAxisLatToShow = axisToShow(4)-axisToShow(3);
    plot3([centerLon centerLon], [axisToShow(3)-deltaAxisLatToShow, ...
        axisToShow(4)+deltaAxisLatToShow], ones(2,1), ...
        '-.w', 'LineWidth', 2);
    
    % Hide Google map for the lidar data.
    if false
        fill3([centerLon centerLon ...
            axisToShow(2)+deltaAxisLatToShow ...
            axisToShow(2)+deltaAxisLatToShow], ...
            [axisToShow(3)-deltaAxisLatToShow ...
            axisToShow(4)+deltaAxisLatToShow ...
            axisToShow(4)+deltaAxisLatToShow ...
            axisToShow(3)-deltaAxisLatToShow], ones(1,4).*0.0001, ...
            'w'); %#ok<UNRCH>
    end
    
    % Foliage area.
    treeHeight = VEG_AREA_IMG_META.ZS - VEG_AREA_IMG_META.ALTS;
    boolsIsVegArea = treeHeight(:)>=VEG_AREA_IMG_META.MIN_TREE_H;
    boolsToShowFoliageAreas ...
        = boolsIsVegArea&(VEG_AREA_IMG_META.LONS(:)<centerLon) ...
        &(VEG_AREA_IMG_META.LATS(:)<centerLat);
    hFoliageAreas = scatter3(...
        VEG_AREA_IMG_META.LONS(boolsToShowFoliageAreas), ...
        VEG_AREA_IMG_META.LATS(boolsToShowFoliageAreas), ...
        ones(sum(boolsToShowFoliageAreas),1).*0.9, 'g.');
    foliageScatterAlpha = 1;
    set(hFoliageAreas, 'MarkerFaceAlpha', foliageScatterAlpha, ...
        'MarkerEdgeAlpha', foliageScatterAlpha);
    
    % Center horizontal divider.
    plot3([axisToShow(1)-deltaAxisLonToShow ...
        axisToShow(2)+deltaAxisLonToShow], ...
        [centerLat centerLat], ones(2,1), ...
        '-.w', 'LineWidth', 2);
    
    % Change to X-Y view.
    view(2);
    legend([hTreeLocs, hFoliageAreas], ...
        'Manually labeled tree trunks', ...
        'Automatically extracted foliage regions', ...
        'Location', 'northwest');
    % Show the figure.
    set(hVegOverview, 'visible','on');

    disp('    Saving plot ...')
    
    % Export a copy for papers.
    pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
        'PostProcessingResults', '1_EpsFigs');
        
    %     saveEpsFigForPaper(hVegOverview, ...
    %         fullfile(pathToSavePaperFigs, '2_0_Overview_For_Veg.eps'));
    
    % We will save a .png file instead because export_fig has trouble
    % exporting this figure as an .eps file.
    set(hVegOverview, 'Color', 'white');    
    export_fig(fullfile(pathToSavePaperFigs, ...
        '2_0_Overview_For_Veg.png'), '-png', '-transparent');
    
    disp('    It is preferred to use export_fig to get an .eps copy by:')
    disp('         saveEpsFigForPaper(hVegOverview, ... ')
    disp('             fullfile(pathToSavePaperFigs, ''2_0_Overview_For_Veg.eps''));')
    disp('    If it does not work, one can still manually adjust the figure and save copies as needed.')
    disp('    After adjusting the figure, please use the command below to save a copy.')
    disp('        hgexport(hVegOverview, fullfile(pathToSavePaperFigs, ...')
    disp('            ''ManualWork'', ''2_0_Overview_For_Veg_Manual.eps''));')
    
    disp('    Done!')
end

% EOF