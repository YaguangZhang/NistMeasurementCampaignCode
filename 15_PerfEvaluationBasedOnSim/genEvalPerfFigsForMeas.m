%GENEVALPERFFIGSFORMEAS Generate figures for comparing the simulation
%results with the measurements.
%
% Yaguang Zhang, Purdue, 10/24/2019

% Create directories if necessary.
if exist(curAbsPathToSavePlots, 'dir')~=7
    mkdir(curAbsPathToSavePlots);
end

% Set what measurement results to keep in the figures.
curContiPathLossesWithGpsInfo = contiPathLossesWithGpsInfo;
[curContiPathLossesExtraInfo.contiTxGains, ...
    curContiPathLossesExtraInfo.contiRxGains] = deal(cell(numOfTracks, 1));
for idxTrack = 1:numOfTracks
    curContiPathLossesWithGpsInfo{idxTrack} ...
        (~boolsToKeepMeas{idxTrack}, :) = [];
    curContiPathLossesExtraInfo ...
        .contiTxGains{idxTrack} = contiPathLossesExtraInfo ...
        .contiTxGains{idxTrack}(boolsToKeepMeas{idxTrack});
    curContiPathLossesExtraInfo ...
        .contiRxGains{idxTrack} = contiPathLossesExtraInfo ...
        .contiRxGains{idxTrack}(boolsToKeepMeas{idxTrack});
end
allContiPathLossesWithGpsInfo = vertcat(curContiPathLossesWithGpsInfo{:});

%% TX Antenna Main Lobe Figs
% For each track's measurement results, generate an overview plot with the
% TX antenna coverage area. We will show both the HPBW and FNBW areas as
% triangle polyshapes.

disp(' ');
disp('Plotting antenna main lobe illustrations ...');

EXPECTED_PATH_LOSS_RANGE = [195, 50]; % [max(allMeas), min(allMeas)]

numOfTicklabels = 11;
figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
    39.9893839683981, 39.9915745444857];

% Overview.
allRxLonLats = allContiPathLossesWithGpsInfo(:, [3,2]);
allMeas = allContiPathLossesWithGpsInfo(:, 1);

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

curDirToSaveFig = fullfile(curAbsPathToSavePlots, ...
    'Overview_BasicPLs.png');
saveas(curFig, curDirToSaveFig);

% By track.
for idxTrack = 1:numOfTracks
    curPLs = curContiPathLossesWithGpsInfo{idxTrack}(:,1);
    curRxLonLats = curContiPathLossesWithGpsInfo{idxTrack}(:,[3,2]);
    curTxAzInDeg = TX_INFO_LOGS{1}(idxTrack).txAz;
    
    
    hpbwPolyshape = hpbwLonLatPolyshapes{idxTrack};
    fnbwPolyshape = fnbwLonLatPolyshapes{idxTrack};
    
    curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
    hTx = plot(TX_LON, TX_LAT, '^g', ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.5);
    colormap hot;
    if ~isempty(curRxLonLats)
        plot3k([curRxLonLats, curPLs], ...
            'Marker', {'.', 12}, ...
            'ColorRange', EXPECTED_PATH_LOSS_RANGE, ...
            'CBLabels', numOfTicklabels, 'Labels', ...
            {['Track #', num2str(idxTrack)], ...
            'Longitude', 'Latitude', ...
            '', 'Path Loss (dB)'});
    end
    xticks([]); yticks([]);
    
    plot(fnbwPolyshape, 'FaceColor', 'cyan');
    plot(hpbwPolyshape, 'FaceColor', 'green');
    
    axis(figAxisToSet); view(2);
    if ~isempty(curRxLonLats)
        plotGoogleMapAfterPlot3k(curFig, 'satellite');
    else
        plot_google_map('MapType', 'satellite');
    end
    
    curDirToSaveFig = fullfile(curAbsPathToSavePlots, ...
        ['basicPLsWithTxMainLobe_Track_', num2str(idxTrack), '.png']);
    saveas(curFig, curDirToSaveFig);
end

%% Load Simulation Results for the Measured Locations

disp(' ');
disp('    Comparing simulation results with measurements...');

% Find all simulation result files.
simXlsxFilesForMeasTracks ...
    = rdir(fullfile(ABS_DIR_TO_SIM_RESULTS, '*.xlsx'));
numOfTracks = length(curContiPathLossesWithGpsInfo);
assert(length(simXlsxFilesForMeasTracks)==numOfTracks, ...
    ['Expecting ', num2str(numOfTracks), ...
    ' simulation result .xlsx files, found ', ...
    num2str(length(simXlsxFilesForMeasTracks)), ' instead!']);

% Verify the filenames and order them according to track numbers.
numOfSimResultsXlsx = length(simXlsxFilesForMeasTracks);
dirsToLoadSimResultsForEachTrack = cell(numOfTracks,1);
for idxXlsx = 1:numOfSimResultsXlsx
    [~, curSimResultFileName] ...
        = fileparts(simXlsxFilesForMeasTracks(idxXlsx).name);
    curSimResultTrackNum = regexp(lower(curSimResultFileName), ...
        'track_(\d+)', 'tokens');
    assert(length(curSimResultTrackNum)==1, ...
        'Only one track index is expected in the filename!');
    dirsToLoadSimResultsForEachTrack{ ...
        str2double(cell2mat(curSimResultTrackNum{1}))} ...
        = simXlsxFilesForMeasTracks(idxXlsx).name;
end

%% For locating outlayers.
numOfSigmasOutlayer = 1;

% Load the simulation results and generate comparison plots.
[curAllShiftedSims, curAllDiffBetweenSimAndMeas] ...
    = deal(cell(numOfTracks,1));
for idxTrack = 1:numOfTracks
    curDirToLoadSimResults = dirsToLoadSimResultsForEachTrack{idxTrack};
    if ~isempty(curDirToLoadSimResults)
        curSimLoss = loadSimLossFromExcel(curDirToLoadSimResults);
        curSimLoss = curSimLoss(boolsToKeepMeas{idxTrack});
        
        assert(all(~isnan(curSimLoss)), ...
            'NaN value found in simulation results!');
        
        curMeasLosses = curContiPathLossesWithGpsInfo{idxTrack}(:,1);
        expectedNumOfSamps = length(curMeasLosses);
        
        [curCalibratedSim, curBestShift, curMultiFactor] ...
            = calibrateSimPlsWithMeas(curSimLoss, curMeasLosses);
        
        % Load history TX to RX distance.
        curRxLocCsv = readtable(fullfile(ABS_PATH_TO_SIM_CSV_FILES, ...
            ['rxLoc_meas_', num2str(idxTrack), '.csv']));
        curRxToTx3DDistInM ...
            = curRxLocCsv.rxToTx3DDistInM(boolsToKeepMeas{idxTrack});
        
        % Plot RMSD vs shift around the best shift value.
        curFigFilenamePrefix = 'SimVsMeas';
        
        bestRmsd = sqrt(mean((curCalibratedSim-curMeasLosses).^2));       
        
        % Plot simulation results with measurements.
        hFigSimVsMeasByIdx = figure('visible', ~flagGenFigSilently);
        hold on;
        hSim = plot(1:expectedNumOfSamps, curCalibratedSim, 'x');
        hMeas = plot(1:expectedNumOfSamps, curMeasLosses, '.');
        title({['shift = ', ...
                num2str(curBestShift, '%.2f'), ' dB, multiFactor = ', ...
                num2str(curMultiFactor, '%.2f')]; ...
                ['Best RMSD = ', num2str(bestRmsd, '%.2f'), ' dB']});
        xlabel('Sample'); ylabel('RMSD (dB)');
        grid on; grid minor; axis tight;
        if ~isempty(hSim)
            legend([hSim, hMeas], 'Simulation', 'Measurement');
        end
        pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
            [curFigFilenamePrefix, '_PlVsSampIdx_Track_', ...
            num2str(idxTrack), '.png']);
        
        saveas(hFigSimVsMeasByIdx, pathToSaveCurFig);
        
        % Another comparison figure with TX-to-RX distance as the x axis.
        [xs, indicesSortByDist] = sort(curRxToTx3DDistInM);
        
        hFigSimVsMeasByDist = figure('visible', ~flagGenFigSilently);
        hold on;
        if ~isempty(xs)
            hSim = plot(xs, curCalibratedSim(indicesSortByDist), 'x-');
            hMeas = plot(xs, curMeasLosses(indicesSortByDist), '.--');
        end
        title({['shift = ', ...
                num2str(curBestShift, '%.2f'), ' dB, multiFactor = ', ...
                num2str(curMultiFactor, '%.2f')]; ...
                ['Best RMSD = ', num2str(bestRmsd, '%.2f'), ' dB']});
        xlabel('3D RX-to-TX Distance (m)'); ylabel('RMSD (dB)');
        grid on; grid minor; axis tight;
        if ~isempty(hSim)
            legend([hSim, hMeas], 'Simulation', 'Measurement', ...
                'Location', 'SouthEast');
        end
        pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
            [curFigFilenamePrefix, '_PlVsDist_Track_', ...
            num2str(idxTrack), '.png']);
        
        saveas(hFigSimVsMeasByDist, pathToSaveCurFig);
        
        % Statistical performance comparison: the RMSD comparing the
        % simulation results with the measurements.
        curDiffBetweenSimAndMeas = curCalibratedSim - curMeasLosses;
        curLats = curContiPathLossesWithGpsInfo{idxTrack}(:, 2);
        curLons = curContiPathLossesWithGpsInfo{idxTrack}(:, 3);
        
        [~, ~, ~, ~, hCurFigFittedNorm] ...
            = fitDataToNormDist(curDiffBetweenSimAndMeas);
        pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
            [curFigFilenamePrefix, '_FittedNormForDiff_', ...
            num2str(idxTrack), '.png']);
        
        saveas(hCurFigFittedNorm, pathToSaveCurFig);
        
        % For showing outlayers on map.
        hCurOutlayersOnMap = plotOutlayersOnMap( ...
            curDiffBetweenSimAndMeas, curLons, curLats, ...
            numOfSigmasOutlayer, figAxisToSet);
        pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
            [curFigFilenamePrefix, '_Outlayers_', ...
            num2str(idxTrack), '.png']);
        
        saveas(hCurOutlayersOnMap, pathToSaveCurFig);
        
        % Store the results for overview plots.
        curAllShiftedSims{idxTrack} = curCalibratedSim;
        curAllDiffBetweenSimAndMeas{idxTrack} = curDiffBetweenSimAndMeas;
    else
        warning(['Missing simulation results for track #', ...
            num2str(idxTrack)]);
    end
end

% The RMSD comparison for all data.
allMeasLosses = allContiPathLossesWithGpsInfo(:,1);
allShiftedSims = vertcat(curAllShiftedSims{:});
[~, ~, ~, ~, hCurFigFittedNorm] ...
    = fitDataToNormDist(allShiftedSims - allMeasLosses);
pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
    [curFigFilenamePrefix, '_FittedNormForDiff_AllMeasLocs.png']);

saveas(hCurFigFittedNorm, pathToSaveCurFig);

% The outlayers for all data.
allDiffBetweenSimAndMeas = vertcat(curAllDiffBetweenSimAndMeas{:});
allMeasLats = allContiPathLossesWithGpsInfo(:, 2);
allMeasLons = allContiPathLossesWithGpsInfo(:, 3);

hCurOutlayersOnMap = plotOutlayersOnMap( ...
    allDiffBetweenSimAndMeas, allMeasLons, allMeasLats, ...
    numOfSigmasOutlayer, figAxisToSet);
pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
    [curFigFilenamePrefix, '_Outlayers_AllMeasLocs.png']);

saveas(hCurOutlayersOnMap, pathToSaveCurFig);

if flagGenFigSilently
    close all;
end

disp('    Done!');
% EOF