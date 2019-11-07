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
curAllContiPathLossesWithGpsInfo ...
    = vertcat(curContiPathLossesWithGpsInfo{:});

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
allContiPathLossesWithGpsInfoOriginal ...
    = vertcat(contiPathLossesWithGpsInfo{:});
allRxLonLats = allContiPathLossesWithGpsInfoOriginal(:, [3,2]);
allMeas = allContiPathLossesWithGpsInfoOriginal(:, 1);
allBoolsToKeepMeas = vertcat(boolsToKeepMeas{:});

curFig = figure('Visible', ~FLAG_GEN_FIGS_SILENTLY); hold on;
hTx = plot(TX_LON, TX_LAT, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
colormap hot;
plot3k([allRxLonLats(allBoolsToKeepMeas, :), ...
    allMeas(allBoolsToKeepMeas)], ...
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

if flagGenFigsForAttennaBeam
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
[curAllShiftedSims, curAllDiffBetweenSimAndMeas, ...
    curAllSimLosses, curAllMeasLosses, curAllRxToTx3DDistInM] ...
    = deal(cell(numOfTracks,1));
for idxTrack = 1:numOfTracks
    curDirToLoadSimResults = dirsToLoadSimResultsForEachTrack{idxTrack};
    if ~isempty(curDirToLoadSimResults)
        curSimLosses = loadSimLossFromExcel(curDirToLoadSimResults, 1);
        curSimLosses = curSimLosses(boolsToKeepMeas{idxTrack});
        
        assert(all(~isnan(curSimLosses)), ...
            'NaN value found in simulation results!');     
        
        curMeasLosses = curContiPathLossesWithGpsInfo{idxTrack}(:,1);
        expectedNumOfSamps = length(curMeasLosses);
        
        % Load history TX to RX distance.
        curRxLocCsv = readtable(fullfile(ABS_PATH_TO_SIM_CSV_FILES, ...
            ['rxLoc_meas_', num2str(idxTrack), '.csv']));
        curRxToTx3DDistInM ...
            = curRxLocCsv.rxToTx3DDistInM(boolsToKeepMeas{idxTrack});
        
        [curSortedRxToTxDists, indicesSortByDist] ...
            = sort(curRxToTx3DDistInM);
        
        curSimLosses = curSimLosses(indicesSortByDist);
        curMeasLosses = curMeasLosses(indicesSortByDist);
        
        [curCalibratedSim, curBestShift, curMultiFactor] ...
            = calibrateSimPlsWithMeas(curSimLosses, curMeasLosses);
        
        curAllSimLosses{idxTrack} = curCalibratedSim;
        curAllMeasLosses{idxTrack} = curMeasLosses;
        curAllRxToTx3DDistInM{idxTrack} = curSortedRxToTxDists;
        
        % Plots.
        curFigFilenamePrefix = 'SimVsMeas';        
        curRmsd = sqrt(mean((curCalibratedSim-curMeasLosses).^2));
        
        if expectedNumOfSamps>0
            % Plot raw simulation results and meassurements in the same
            % plot.
            hFigRawSimVsMeasByDist ...
                = figure('visible', ~flagGenFigSilently);
            hold on;
            yyaxis left;
            hSim = plot(curSortedRxToTxDists, curSimLosses, 'x');
            xlabel('3D RX-to-TX Distance (m)'); 
            ylabel('Sim Peak Diff (dB)');
            yyaxis right;
            hMeas = plot(curSortedRxToTxDists, curMeasLosses, '.');
            ylabel('Path Loss (dB)');            
            grid on; grid minor; 
            pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
                [curFigFilenamePrefix, '_RawSimVsMeasByDist_Track_', ...
                num2str(idxTrack), '.png']);
            
            saveas(hFigRawSimVsMeasByDist, pathToSaveCurFig);
            
            % Plot simulation results with measurements.
            hFigSimVsMeasByIdx = figure('visible', ~flagGenFigSilently);
            hold on;
            hSim = plot(1:expectedNumOfSamps, curCalibratedSim, 'x');
            hMeas = plot(1:expectedNumOfSamps, curMeasLosses, '.');
            title({['shift = ', ...
                num2str(curBestShift, '%.2f'), ' dB, multiFactor = ', ...
                num2str(curMultiFactor, '%.2f')]; ...
                ['Best RMSD = ', num2str(curRmsd, '%.2f'), ' dB']});
            xlabel('Sample'); ylabel('Path Loss (dB)');
            grid on; grid minor; axis tight;
            if ~isempty(hSim)
                legend([hSim, hMeas], 'Simulation', 'Measurement');
            end
            pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
                [curFigFilenamePrefix, '_PlVsSampIdx_Track_', ...
                num2str(idxTrack), '.png']);
            
            saveas(hFigSimVsMeasByIdx, pathToSaveCurFig);
        end
        
        if ~isempty(curSortedRxToTxDists)
            % Another comparison figure with TX-to-RX distance as the x
            % axis.
            hFigSimVsMeasByDist = figure('visible', ~flagGenFigSilently);
            hold on;
            hSim = plot(curSortedRxToTxDists, ...
                curCalibratedSim, 'x-');
            hMeas = plot(curSortedRxToTxDists, ...
                curMeasLosses, '.--');
            title({['shift = ', ...
                num2str(curBestShift, '%.2f'), ' dB, multiFactor = ', ...
                num2str(curMultiFactor, '%.2f')]; ...
                ['Best RMSD = ', num2str(curRmsd, '%.2f'), ' dB']});
            xlabel('3D RX-to-TX Distance (m)'); ylabel('Path Loss (dB)');
            grid on; grid minor; axis tight;
            if ~isempty(hSim)
                legend([hSim, hMeas], 'Simulation', 'Measurement', ...
                    'Location', 'SouthEast');
            end
            pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
                [curFigFilenamePrefix, '_PlVsDist_Track_', ...
                num2str(idxTrack), '.png']);
            
            saveas(hFigSimVsMeasByDist, pathToSaveCurFig);
        end
        
        % Statistical performance comparison: the RMSD comparing the
        % simulation results with the measurements.
        curDiffBetweenSimAndMeas = curCalibratedSim - curMeasLosses;
        if ~isempty(curDiffBetweenSimAndMeas)
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
        end
        
        % Store the results for overview plots.
        curAllShiftedSims{idxTrack} = curCalibratedSim;
        curAllDiffBetweenSimAndMeas{idxTrack} = curDiffBetweenSimAndMeas;
    else
        warning(['Missing simulation results for track #', ...
            num2str(idxTrack)]);
    end
end

curSimLosses = vertcat(curAllSimLosses{:});
curMeasLosses = vertcat(curAllMeasLosses{:});
curRxToTx3DDistInM = vertcat(curAllRxToTx3DDistInM{:});

[curSortedRxToTxDists, indicesSortByDist] = sort(curRxToTx3DDistInM);
        
curCalibratedSim = curSimLosses(indicesSortByDist);
curMeasLosses = curMeasLosses(indicesSortByDist);

curRmsd = sqrt(mean((curCalibratedSim-curMeasLosses).^2));

% The path loss over distance plot for all data.
hFigSimVsMeasByDist = figure('visible', ~flagGenFigSilently);
hold on;
hSim = plot(curSortedRxToTxDists, ...
    curCalibratedSim, 'x-');
hMeas = plot(curSortedRxToTxDists, ...
    curMeasLosses, '.--');
title(['Overall RMSD = ', num2str(curRmsd, '%.2f'), ' dB']);
xlabel('3D RX-to-TX Distance (m)'); ylabel('Path Loss (dB)');
grid on; grid minor; axis tight;
if ~isempty(hSim)
    legend([hSim, hMeas], 'Simulation', 'Measurement', ...
        'Location', 'SouthEast');
end
pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
    [curFigFilenamePrefix, '_PlVsDist_AllTracks.png']);

saveas(hFigSimVsMeasByDist, pathToSaveCurFig);
        
% The RMSD comparison for all data.
allMeasLosses = curAllContiPathLossesWithGpsInfo(:,1);
allShiftedSims = vertcat(curAllShiftedSims{:});
[~, ~, ~, ~, hCurFigFittedNorm] ...
    = fitDataToNormDist(allShiftedSims - allMeasLosses);
pathToSaveCurFig = fullfile(curAbsPathToSavePlots, ...
    [curFigFilenamePrefix, '_FittedNormForDiff_AllMeasLocs.png']);

saveas(hCurFigFittedNorm, pathToSaveCurFig);

% The outlayers for all data.
allDiffBetweenSimAndMeas = vertcat(curAllDiffBetweenSimAndMeas{:});
allMeasLats = curAllContiPathLossesWithGpsInfo(:, 2);
allMeasLons = curAllContiPathLossesWithGpsInfo(:, 3);

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