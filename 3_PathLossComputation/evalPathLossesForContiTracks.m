% EVALPATHLOSSESFORCONTITRACKS Load the Conti track measurements and their
% cooresponding GPS logs form the NIST dataset, and generate for each track
% a path loss over track (on map) plot.
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

% We will need the function genCalibrationFct.m for calibration.
addpath(fullfile(pwd, '1_Calibration'));

% The functions:
%   -  thresholdWaveform.m for noise elimination
%    - antPatInter.m for antenna gain calculation
% have been included in lib.

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputationConti');
% There is one case where we recorded the reflection of a spot without
% moving the RX at all at USNA. The Tx info log for that one is left empty.
ABS_PATH_TO_SAVE_PLOTS_REF = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'reflection');

% Reuse results from plotInfo.m, calibrateRx.m, fetchAntennaPattern.m, and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_DATA_INFO_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'SummaryReport', 'plots', 'plotInfo.mat');
ABS_PATH_TO_CALI_LINES_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration', 'lsLinesPolys.mat');
% We have manually copied antennaPatternFlipped.mat from the EARS dataset.
ABS_PATH_TO_ANT_PAT_FILE = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'AntennaPattern', ...
    'antennaPatternFlipped.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

% For setting the threshold during the noise elimination.
% NUM_SIGMA_FOR_THRESHOLD = 3.5;

% We will use the Google service for RX altitudes.
FLAG_USE_GOOGLE_FOR_ALT = true;
GOOGLE_MAPS_API = 'AIzaSyDlkaE_QJxvRJpTutWG0N-LCvoT0e7FPHE';

% Set this flag to be true to skip the computing process and reuse history
% results for generating plots only.
FLAG_ONLY_GEN_FIGS = false;

% Manually set the figure size and axis to best show the data on plots.
figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
    39.9893839683981, 39.9915745444857];

%% Before Processing the Data

disp(' -------------------------- ')
disp('  pathLossesForContiTracks')
disp(' -------------------------- ')

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS_REF, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS_REF);
end

%% Get Info for Measurement Data Files and Calibration Polynomials

disp(' ')
disp('    Loading results from: ')
disp('      - plotInfo.m')
disp('      - calibrateRx.m')
disp('      - antennaPattern.m')
disp('      - loadMeasCampaignInfo.m')

assert(exist(ABS_PATH_TO_DATA_INFO_FILE, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');
assert(exist(ABS_PATH_TO_CALI_LINES_FILE, 'file')==2, ...
    'Couldn''t find lsLinesPolys.mat! Please run PostProcessing/1_Calibration/calibrateRx.m first.');
assert(exist(ABS_PATH_TO_ANT_PAT_FILE, 'file')==2, ...
    'Couldn''t find antennaPatternFlipped.mat! Please manually copy it to PostProcessingResults/AntennaPattern/.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/3_PathLossComputation/loadMeasCampaignInfo.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'allSeriesParentDirs' and 'allSeriesDirs'.
load(ABS_PATH_TO_DATA_INFO_FILE);
% Get 'lsLinesPolys', 'lsLinesPolysInv', 'fittedMeaPs', 'fittedCalPs', and
% 'rxGains'.
load(ABS_PATH_TO_CALI_LINES_FILE);
% Get 'pat28GAzNorm', and 'pat28GElNorm'.
load(ABS_PATH_TO_ANT_PAT_FILE);

% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Sample rate used for GnuRadio. Needed by computePathLossForOutFileDir.m.
Fs = F_S;
% TX power (after upconverter) in dBm.
txPower  = TX_POWER_DBM;

disp('    Done!')

%% Search for the Conti Measurement Data Files

if FLAG_ONLY_GEN_FIGS
    disp(' ')
    disp('    FLAG_ONLY_GEN_FIGS is set to be true ...')
    disp('    We will load history results and only generate figures ...')
    disp(' ')
else
    disp(' ')
    disp('    Searching for measurement data files ...')
    
    % Here we don't care about when the data were collected so let's
    % rearrange all the dir struct into one array.
    allSeriesDirsArray = vertcat(allSeriesDirs{1:end});
    numSeries = length(allSeriesDirsArray);
    
    % Scan the series folder for Gnu Radio out files, as well as the
    % corresponding GPS log files.
    [contiOutFilesDirs, contiGpsFilesDirs] = deal(cell(numSeries,1));
    
    for idxContiOutFile = 1:numSeries
        disp(['        Scanning series folder ', num2str(idxContiOutFile), ...
            '/', num2str(numSeries), '...']);
        
        % Here we only load the _Conti .out files.
        regexPattern = '\d+_NistFoliage';
        [contiOutFilesDirs{idxContiOutFile}, contiGpsFilesDirs{idxContiOutFile}] = ...
            loadFilesFromSeriesDir(allSeriesDirsArray(idxContiOutFile), regexPattern);
    end
    
    % Get rid of empty elements.
    boolsEmptyDirs = cellfun(@(d) isempty(d), contiOutFilesDirs);
    contiOutFilesDirs(boolsEmptyDirs) = [];
    contiGpsFilesDirs(boolsEmptyDirs) = [];
    
    disp('    Done!')
end

%% Compute Path Losses

pathContiPathLossesFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'contiPathLossesWithGpsInfo.mat');

if FLAG_ONLY_GEN_FIGS
    disp('    Loading history results ...')
    if ~exist(pathContiPathLossesFileToSave, 'file')
        error(['History path loss results is not available! ', ...
            'Please run this script with FLAG_ONLY_GEN_FIGS=false to generate them.'])
    end
    load(pathContiPathLossesFileToSave);
else
    disp(' ')
    disp('    Computing path losses...')
    
    % Compute the path losses and save them into a matrix together with the
    % GPS info.
    numContiOutFiles = sum(cellfun(@(d) length(d), contiOutFilesDirs));
    % More specifically, we will create one cell for each conti track to
    % store a matrix, where each row is a [path loss (dB), lat, lon, alt]
    % array for a GPS position.
    [contiPathLossesWithGpsInfo, contiPathLossesExtraInfoTxGains, ...
        contiPathLossesExtraInfoRxGains, contiPathLossesExtraInfoOriAlts, ...
        contiPathLossesExtraInfoSampIndices, contiPathLossesExtraLocInfo] ...
        = deal(cell(numContiOutFiles, 1));
    % Also save the meta info needed to map the path loss back to the
    % measurements. We choose to save the full file path to the .out file
    % for convenience.
    absPathsContiOutFiles = cell(numContiOutFiles, 1);
    % For indexing and storing the results into cells.
    contiOutFileCounter = 1;
    contiOutFileIndicesReflection = [];
    
    for idxContiOutFile = 1:numContiOutFiles
        curContiOutFileDir = contiOutFilesDirs{idxContiOutFile};
        curContiGpsFilesDirs = contiGpsFilesDirs{idxContiOutFile};
        
        disp('-------------------------------------------------')
        disp(['        Processing conti. out file ', ...
            num2str(idxContiOutFile), '/', ...
            num2str(numContiOutFiles), '...']);
        disp('-------------------------------------------------')
        disp(['            Full path: ', ...
            fullfile(curContiOutFileDir.folder, curContiOutFileDir.name)]);
        
        % Compute path losses (without considering the antenna gain) for
        % this track. We will use the amplitude version of
        % thresholdWaveform.m without plots for debugging as the noise
        % eliminiation function.
        noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
        [ curContiPathLossesWithGpsInfo, curAbsPathOutFile, ...
            curContiOriAlts, curOutFileSampleRanges ] ...
            = computePathLossesForContiOutFileDir(...
            curContiOutFileDir, curContiGpsFilesDirs, ...
            noiseEliminationFct, FLAG_USE_GOOGLE_FOR_ALT, GOOGLE_MAPS_API);
        
        % Fetch the measurement campaign meta records.
        [absCurParDir, curSeries] = fileparts(curContiOutFileDir.folder);
        idxParDir = find(cellfun(@(d) strcmp(d, absCurParDir), ...
            TX_INFO_LOGS_ABS_PAR_DIRS));
        idxCurSer = str2double(curSeries((length('Series_')+1):end));
        assert(length(idxParDir)==1, ...
            'Error: More than 1 parent folder match with the meta information of the current Series data!');
        curTxInfoLog = TX_INFO_LOGS{idxParDir}(idxCurSer);
        assert(curTxInfoLog.seriesNum==idxCurSer, ...
            'The series index number in the matched Tx info log does not agree with idxCurSer.');
        
        [numPathLosses, ~] = size(curContiPathLossesWithGpsInfo);
        [curContiTxGains, curContiRxGains] ...
            = deal(nan(numPathLosses, 1));
        curContiExtraLocInfo = nan(numPathLosses, 7);
        % For the cases where TX orientation info is nan, we will ignore
        % the antenna gain calculation.
        if ~any(isnan([curTxInfoLog.txAz, curTxInfoLog.txEl, ...
                curTxInfoLog.rxAz, curTxInfoLog.rxEl]))
            disp(' ')
            disp('            Computing antenna gains in 3D...')
            for idxPathLoss = 1:numPathLosses
                % Compute the antenna gains.
                %    function [ gain ] = computeAntGain(lat0, lon0, alt0,
                %    ...
                %        az0, el0, ...
                %         lat, lon, alt, ...
                %        antPatAz, antPatEl)
                txGain = computeAntGain(TX_LAT, TX_LON, TX_ALT+TX_HEIGHT_M, ...
                    curTxInfoLog.txAz, curTxInfoLog.txEl, ...
                    curContiPathLossesWithGpsInfo(idxPathLoss,2), ... % Rx lat
                    curContiPathLossesWithGpsInfo(idxPathLoss,3), ... % Rx lon
                    curContiPathLossesWithGpsInfo(idxPathLoss,4) ... % Rx alt
                    +RX_HEIGHT_M, ...
                    pat28GAzNorm, pat28GElNorm);
                rxGain = computeAntGain( ...
                    curContiPathLossesWithGpsInfo(idxPathLoss,2), ... % Rx lat
                    curContiPathLossesWithGpsInfo(idxPathLoss,3), ... % Rx lon
                    curContiPathLossesWithGpsInfo(idxPathLoss,4) ... % Rx alt
                    +RX_HEIGHT_M, ...
                    curTxInfoLog.rxAz, curTxInfoLog.rxEl, ...
                    TX_LAT, TX_LON, TX_ALT+TX_HEIGHT_M, ...
                    pat28GAzNorm, pat28GElNorm);
                
                disp(['            Sample ', num2str(idxPathLoss), ...
                    '/', num2str(numPathLosses), ...
                    ': txGain = ', num2str(txGain), ...
                    ' dB; rxGain = ', num2str(rxGain), ' dB'])
                
                finalPathLossValue ...
                    = curContiPathLossesWithGpsInfo(idxPathLoss,1) ...
                    + txGain + rxGain;
                
                disp(['                PathLossBefore = ', ...
                    num2str(curContiPathLossesWithGpsInfo(idxPathLoss,1)), ...
                    ' dB; FinalPathLoss = ', ...
                    num2str(finalPathLossValue), ' dB'])
                
                curContiPathLossesWithGpsInfo(idxPathLoss,1) ...
                    = finalPathLossValue;
                curContiTxGains(idxPathLoss) = txGain;
                curContiRxGains(idxPathLoss) = rxGain;
                curContiExtraLocInfo(idxPathLoss, :) ...
                    = [curTxInfoLog.rxAz, curTxInfoLog.rxEl, ...
                    TX_LAT, TX_LON, TX_HEIGHT_M, ...
                    curTxInfoLog.txAz, curTxInfoLog.txEl];
            end
        else
            warning('We not have enough information to extract antenna gains!');
            
            contiOutFileIndicesReflection = [contiOutFileIndicesReflection, ...
                idxContiOutFile]; %#ok<AGROW>
        end
        % Store the results.
        contiPathLossesWithGpsInfo{idxContiOutFile} ...
            = curContiPathLossesWithGpsInfo;
        contiPathLossesExtraInfoTxGains{idxContiOutFile} = curContiTxGains;
        contiPathLossesExtraInfoRxGains{idxContiOutFile} = curContiRxGains;
        contiPathLossesExtraInfoOriAlts{idxContiOutFile} = curContiOriAlts;
        contiPathLossesExtraInfoSampIndices{idxContiOutFile} ...
            = curOutFileSampleRanges;
        contiPathLossesExtraLocInfo{idxContiOutFile} ...
            = curContiExtraLocInfo;
        absPathsContiOutFiles{idxContiOutFile} = curAbsPathOutFile;
    end
    assert(all(~isempty(vertcat(contiPathLossesWithGpsInfo{1:end}))));
    
    disp('    Saving the results...')
    % For absPathsOutFiles, convert it to relative paths under the data
    % folder, which will already contain enough information we need.
    contiOutFilesRelPathsUnderDataFolder = ...
        cellfun(@(p) regexp(p, 'Data[\/\\]([a-zA-Z\d\/\\_]+.out)$', ...
        'tokens'), absPathsContiOutFiles);
    % Extra info for debugging.
    contiPathLossesExtraInfo.contiTxGains = contiPathLossesExtraInfoTxGains;
    contiPathLossesExtraInfo.contiRxGains = contiPathLossesExtraInfoRxGains;
    contiPathLossesExtraInfo.FLAG_USE_GOOGLE_FOR_ALT = FLAG_USE_GOOGLE_FOR_ALT;
    contiPathLossesExtraInfo.contiOriAlts ...
        = contiPathLossesExtraInfoOriAlts;
    contiPathLossesExtraInfo.contiSampIndices ...
        = contiPathLossesExtraInfoSampIndices;
    contiPathLossesExtraInfo.contiExtraLocInfo ...
        = contiPathLossesExtraLocInfo;
    save(pathContiPathLossesFileToSave, ...
        'contiPathLossesWithGpsInfo', ...
        'contiPathLossesExtraInfo', ...
        'contiOutFilesRelPathsUnderDataFolder', ...
        'contiOutFileIndicesReflection');
end

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

numTracks = length(contiPathLossesWithGpsInfo);
for idxTrack = 1:numTracks
    
    curRelPath = cell2mat(contiOutFilesRelPathsUnderDataFolder{idxTrack});
    disp(['        Track ', num2str(idxTrack), '/', ...
        num2str(numTracks), ' - ', ...
        curRelPath, '...']);
    
    % Get the date, series number and timestamp for this track, and use
    % them to construct the filename.
    curFileName = strrep(strrep( ...
        strrep(strrep( curRelPath, '\', '_' ),'/','_'), ...
        'measureSignalCont','Epoch'), '.out', '');
    
    curPathLossesWithGpsInfo = contiPathLossesWithGpsInfo{idxTrack};
    
    boolsInvalidCoor = curPathLossesWithGpsInfo(:,2)==0 ...
        & curPathLossesWithGpsInfo(:,3)==0;
    if any(boolsInvalidCoor)
        warning([num2str(sum(boolsInvalidCoor)), ...
            ' invalid (lat, lon) pairs detected (both are 0).', ...
            ' We will ignore these points together with their path losses.']);
    end
    pathLossesWithValidGps = curPathLossesWithGpsInfo(~boolsInvalidCoor,:);
    
    boolsInfPathloss = isinf(pathLossesWithValidGps(:,1));
    if any(boolsInfPathloss)
        warning([num2str(sum(boolsInfPathloss)), ...
            ' inf path loss detected.', ...
            ' We will show these points at z=0 with different markers.']);
    end
    validPathLossesWithValidGps = pathLossesWithValidGps(~boolsInfPathloss,:);
    
    % Plot path losses on map.
    hPathLossesOnMap = figure; hold on; colormap jet;
    plot(validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), 'w.');
    plot(pathLossesWithValidGps(boolsInfPathloss,3), ...
        pathLossesWithValidGps(boolsInfPathloss,2), 'kx');
    hTx = plot3(TX_LON, TX_LAT, TX_HEIGHT_M, '^w', 'MarkerFaceColor', 'b');
    plot_google_map('MapType','satellite');
    plot3k([validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), ...
        validPathLossesWithValidGps(:,1)], 'Marker', {'.', 12});
    % The command plot_google_map messes up the color legend of plot3k, so
    % we will have to fix it here.
    hCb = findall( allchild(hPathLossesOnMap), 'type', 'colorbar');
    hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));
    hold off; grid on; view(0, 90); legend(hTx, 'TX');
    curTitleLabel = strrep(curFileName, '_', '-');
    if ismember(idxTrack, contiOutFileIndicesReflection)
        title({'Path Losses on Map (Reflection)', curTitleLabel});
    else
        title({'Path Losses on Map (Conti. Track)', curTitleLabel});
    end
    xlabel('Lon'); ylabel('Lat'); zlabel('Path Loss (dB)');
    
    % Plot path losses over distance from Tx.
    validPLWithValidGPSCell = num2cell(validPathLossesWithValidGps, 2);
    distsFromTx = cellfun(@(s) ...
        norm([1000.*lldistkm([s(2) s(3)],[TX_LAT,TX_LON]), ...
        TX_ALT+TX_HEIGHT_M-s(4)]), ...
        validPLWithValidGPSCell);
    
    hPathLossesOverDist = figure; hold on; colormap jet;
    plot3k([distsFromTx, zeros(length(distsFromTx),1), ...
        validPathLossesWithValidGps(:,1)], 'Marker', {'.', 6});
    curAxis = axis;
    axis([min([distsFromTx; 1]), max(distsFromTx)+100, curAxis(3:6)]);
    view(0, 0); set(gca, 'XScale', 'log'); grid on;
    newXTicks = [1,10,100,200,500,1000];
    set(gca, 'XTickLabels',newXTicks);
    set(gca, 'XTick',newXTicks);
    title('Path Losses over Distance (Conti)');
    xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
    
    % Save the plots.
    if ismember(idxTrack, contiOutFileIndicesReflection)
        pathPathossesOnMapFileToSave = fullfile(...
            ABS_PATH_TO_SAVE_PLOTS_REF, ...
            curFileName);
    else
        pathPathossesOnMapFileToSave = fullfile(...
            ABS_PATH_TO_SAVE_PLOTS, ...
            curFileName);
    end
    saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '_map.png']);
    saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '_map.fig']);
    saveas(hPathLossesOverDist, [pathPathossesOnMapFileToSave, '_dist.png']);
    saveas(hPathLossesOverDist, [pathPathossesOnMapFileToSave, '_dist.fig']);
end

disp('    Done!')
disp('')

%% Overview
% Plot all non-reflection tracks on the same plot.
disp('    Plotting all non-reflection tracks on the same map...');

hNonRefPathLossesOnMap = figure; hold on; colormap jet;
hTx = plot3(TX_LON, TX_LAT, TX_HEIGHT_M, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
% Cache all the non-reflection path loss records for plotting.
validNonRefPathLossesWithValidGps = [];
for idxTrack = 1:numTracks
    if ismember(idxTrack, contiOutFileIndicesReflection)
        disp(['        Skipping reflection track ', num2str(idxTrack)])
    else
        curPathLossesWithGpsInfo = contiPathLossesWithGpsInfo{idxTrack};
        
        [boolsValidPathlosses, ...
            boolsInvalidCoor, boolsInfPathloss, boolsInvalidData]  ...
            = checkValidityOfPathLossesWithGpsInfo(curPathLossesWithGpsInfo);
        
        validPathLossesWithValidGps ...
            = curPathLossesWithGpsInfo(boolsValidPathlosses,:);
        infPathLossesWithValidGps = curPathLossesWithGpsInfo( ...
            (~boolsInvalidData) & (~boolsInvalidCoor) & boolsInfPathloss,:);
        
        validNonRefPathLossesWithValidGps ...
            = [validNonRefPathLossesWithValidGps;...
            validPathLossesWithValidGps];
        
        plot(infPathLossesWithValidGps(:,3), ...
            infPathLossesWithValidGps(:,2), 'w.');
        plot(infPathLossesWithValidGps(:,3), ...
            infPathLossesWithValidGps(:,2), 'kx');
    end
end

% We need a smaller figure.
curPos = get(hNonRefPathLossesOnMap, 'Position');
overviewFigSize = 'normal'; % 'small'
if strcmp(overviewFigSize, 'small')
    set(hNonRefPathLossesOnMap, 'Position', [curPos(1:2), 330 260]);
end

% Make sure the color map is intuitive and works in grey scale.
colormap hot;
numOfTicklabels = 11;
if strcmp(overviewFigSize, 'small')
    numOfTicklabels = 7;
end
[~, ~, hPlot3kColorBar] ...
    = plot3k([validNonRefPathLossesWithValidGps(:,3), ...
    validNonRefPathLossesWithValidGps(:,2), ...
    validNonRefPathLossesWithValidGps(:,1)], 'Marker', {'.', 12}, ...
    'ColorRange', [max(validNonRefPathLossesWithValidGps(:,1)) ...
    min(validNonRefPathLossesWithValidGps(:,1))], ...
    'CBLabels', numOfTicklabels, 'Labels', ...
    {'', 'Longitude', 'Latitude', '', ''});
xticks([]); yticks([]);

% Draw a square to show the area where we are gonna illustrate the manually
% labeled trees.
axisValuesForTreeOverview = [-105.27701930, -105.27483020, ...
    39.99008676, 39.99156435];
plot3(axisValuesForTreeOverview([1 1 2 2 1]), ...
    axisValuesForTreeOverview([3 4 4 3 3]), ...
    ones(1,5).*max(validNonRefPathLossesWithValidGps(:,1))+1, ...
    '-.', 'LineWidth', 1, 'Color', ones(3,1));

% Manually adjust the visible area.
axis(figAxisToSet);
plot_google_map('MapType','satellite');

% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hNonRefPathLossesOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

hold off; view(2);
legend(hTx, 'TX', 'Location','southeast');
title('Basic Transmission Losses (dB) on Map');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Path Loss (dB)');

% Repeat for a better output result: Manually adjust the visible area.
axis(figAxisToSet);
plot_google_map('MapType','satellite');

% Save the plot.
fileNamePathossesOnMap = 'allContiTracks_normalSize';
pathPathossesOnMapFileToSave = fullfile(...
    ABS_PATH_TO_SAVE_PLOTS, ...
    fileNamePathossesOnMap);
saveas(hNonRefPathLossesOnMap, [pathPathossesOnMapFileToSave, '_map.png']);
saveas(hNonRefPathLossesOnMap, [pathPathossesOnMapFileToSave, '_map.fig']);

% Export an .eps copy for papers.
pathToSavePaperFigs = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', '1_EpsFigs');
saveEpsFigForPaper(hNonRefPathLossesOnMap, ...
    fullfile(pathToSavePaperFigs, ['3_1_', fileNamePathossesOnMap, '.eps']));

disp('    Done!')

% EOF