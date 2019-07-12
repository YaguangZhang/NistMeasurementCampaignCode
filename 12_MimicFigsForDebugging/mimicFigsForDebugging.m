% MIMICFIGSFORDEBUGGING Generate a data set and plot figures for debugging.
%
% Please try running this script first and if it fails, one needs to get a
% Google Maps API key from
%   https://developers.google.com/maps/documentation/javascript/get-api-key
% and run (after the libs are loaded into the current
% workspace):
%     plot_google_map('APIKey', 'xxxxxxxxxx');
% where xxxxxxxxxx is the API key. Then try running this script again.
%
% Yaguang Zhang, Purdue, 07/11/2019

% Add libs to current path.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(genpath(fullfile(pwd, 'RequiredLibCopies')));

% We will save the RX data into a matrix called dataForSimulationFigs with
% columns being:
%   [lat, lon, pathloss, numOfTrees, distTxToRxIn3D, FSPL].
DIR_TO_TX_DATA = fullfile(pwd, 'txData.mat');
DIR_TO_RX_DATA = fullfile(pwd, 'rxDataMatrix.mat');

%% Generate/Load Data

FLAG_DATA_AVAILABLE = true;
if ~FLAG_DATA_AVAILABLE
    % For Tx.    
    save(DIR_TO_TX_DATA, 'TX_LAT', 'TX_LON');
    
    % For Rx.
    contiPathLossesWithGpsInfoMat = vertcat(contiPathLossesWithGpsInfo{:});
    numsOfTreesInFirstFresnelMat ...
        = vertcat(numsOfTreesInFirstFresnel{:});
    txRxLosDistsMat = vertcat(txRxLosDists{:});
    freeSpacePathLossesMat = vertcat(freeSpacePathLosses{:});
    
    % It is required to manually import all required variables first.
    dataForSimulationFigs ...
        = [contiPathLossesWithGpsInfoMat(:, 2:3), ...
        contiPathLossesWithGpsInfoMat(:,1) numsOfTreesInFirstFresnelMat ...
        txRxLosDistsMat freeSpacePathLossesMat];  
    
    save(DIR_TO_RX_DATA, 'dataForSimulationFigs');
else
    load(DIR_TO_TX_DATA);
    load(DIR_TO_RX_DATA);
end

% Extract the data for readability.
rxLats = dataForSimulationFigs(:, 1);
rxLons = dataForSimulationFigs(:, 2);
rxPathlosses = dataForSimulationFigs(:, 3);
rxNumOfTrees = dataForSimulationFigs(:, 4);
rxDistsFromTxIn3D = dataForSimulationFigs(:, 5);
rxFSPLs = dataForSimulationFigs(:, 6);

%% Fig 1: All Path Loss Values on Map

hPathLossesOnMap = figure; hold on;
hTx = plot(TX_LON, TX_LAT, '^g', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);

% Make sure the color map is intuitive and works in grey scale.
colormap hot;
numOfTicklabels = 11;
[~, ~, hPlot3kColorBar] ...
    = plot3k([rxLons, rxLats, rxPathlosses], 'Marker', {'.', 12}, ...
    'ColorRange', [max(rxPathlosses) ...
    min(rxPathlosses)], ...
    'CBLabels', numOfTicklabels, 'Labels', ...
    {'Basic Transmission Losses (dB) on Map', 'Longitude', 'Latitude', ...
    '', 'Path Loss (dB)'});
xticks([]); yticks([]);

% Manually adjust the visible area.
figAxisToSet = [-105.2774429259207, -105.2744429246357, ...
    39.9893839683981, 39.9915745444857];
axis(figAxisToSet);
plot_google_map('MapType','satellite');

% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hPathLossesOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

hold off; view(2);
legend(hTx, 'TX', 'Location','southeast');

% Repeat for a better output result: Manually adjust the visible area.
axis(figAxisToSet);
plot_google_map('MapType','satellite');

%% Fig 2: Path Loss vs Dist Colored with Tree Numbers

hMeasVsLogDistSimpColoredByTreeNum ...
    = figure('Unit', 'pixels');
% Enlarge the figure horizontally a little bit.
curFigPos = get(hMeasVsLogDistSimpColoredByTreeNum, ...
    'Position');
curFigPos(3) = 640;
curFigPos(4) = 420;
set(hMeasVsLogDistSimpColoredByTreeNum, ...
    'Position', curFigPos);

hold on;
% Measured path losses.
measMarkerSize = 3;
[sortedLosDist, indicesForSortedDist] = sort(rxDistsFromTxIn3D);
hMeas = plot3(sortedLosDist, ...
    rxPathlosses(indicesForSortedDist) - rxFSPLs(indicesForSortedDist), ...
    rxNumOfTrees(indicesForSortedDist), ...
    'ok', 'MarkerSize', measMarkerSize);
maxRxNumOfTrees = max(rxNumOfTrees);
[~, hPlot3kAxis, hPlot3kColorbar] = plot3k([sortedLosDist, ...
    rxPathlosses(indicesForSortedDist) - rxFSPLs(indicesForSortedDist), ...
    rxNumOfTrees(indicesForSortedDist)], ...
    'Marker', {'o', measMarkerSize}, 'Labels', ...
    {'', 'Distance to TX (m)', 'Excess Path Loss (dB)', '', ''}, ...
    'CBFormat', '%2.0f', 'CBLabels', maxRxNumOfTrees+1);
% Plot3k does not label the colorbar correctly. Fix it manually.
hPlot3kColorbar.Ticks ...
    = (0:maxRxNumOfTrees)./maxRxNumOfTrees;
% Verticle colorbar title.
ylabel(hPlot3kColorbar, 'Number of Trees in the 1st Fresnel Zone ');
% Remove marker face color.
hPlot3kLines = findall(hPlot3kAxis, 'type', 'line');
arrayfun(@(hLine) set(hLine, 'MarkerFaceColor', 'none'), hPlot3kLines);

hLegend = legend(hMeas, ...
    'Measurements colored by tree number', ...
    'Location', 'northwest');

axis tight; transparentizeCurLegends; view(2); set(gca, 'xscale','log');
grid on; grid minor;
% Manually set the x tick labels.
xticks([10, 20, 30, 50, 100, 200]);
xticklabels(arrayfun(@(x) {num2str(x)}, xticks'));

% Manually move the legend a little bit to avoid blocking plots.
set(hLegend, 'Units', 'pixel');
curLegendPos = get(hLegend, 'Position');
legPosToSet = [curLegendPos(1)-5, curLegendPos(2), curLegendPos(3:4)];
set(hLegend, 'Position', legPosToSet);

% EOF