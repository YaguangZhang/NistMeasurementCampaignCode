function [ hMapFig, indicesRelavantGpsSamps ] = locateSigSegOnMap(...
    contiPathLossesWithGpsInfo, contiPathLossesExtraInfo, ...
    idxTrack, sampIndexRange, fullPathToSavePlot, axisToSet)
%LOCATESIGSEGONMAP Locate the signal segment in a continuous track on a
%map.
%
%   Inputs:
%   	- contiPathLossesWithGpsInfo, contiPathLossesExtraInfo
%         Path loss computation restults from
%         3_PathLossComputation/evalPathLossesForContiTracks.m.
%       - idxTrack
%         The index for the track from which the signal segment of interest
%         has been extracted.
%       - sampIndexRange
%         The sample index range in the track for the signal segment of
%         interest.
%       - fullPathToSavePlot
%         Optional. The full path to save the plot. If present, the figure
%         will be saved into file and then closed automatically.
%       - axisToSet
%         Optional. If present, the visible area on the map plot will be
%         set accordingly.
%
%   Outputs:
%       - hMapFig
%         The map plot generated to show the singal segment's location.
%       - indicesRelavantGpsSamps
%         The indices for the corresponding GPS track to locate where the
%         plotted signal was recorded.
%
% Yaguang Zhang, Purdue, 06/28/2018

[numSegs, ~] = size(contiPathLossesExtraInfo.contiSampIndices{idxTrack});
indicesRelavantGpsSamps = ...
    find(arrayfun(@(idx) ...
    ~isempty(range_intersection(...
    contiPathLossesExtraInfo.contiSampIndices{idxTrack}(idx,:), ...
    sampIndexRange)), ...
    1:numSegs));

contiPathLossesWithGpsInfoMatrix = vertcat(contiPathLossesWithGpsInfo{:});
contiPathLossesWithGpsInfoFirstPts ...
    = cell2mat(cellfun(@(m) m(1,:), contiPathLossesWithGpsInfo, ...
    'UniformOutput', false));

hMapFig = figure; hold on;
% Plot all the GPS samples.
hAllGps = plot(contiPathLossesWithGpsInfoMatrix(:,3), ...
    contiPathLossesWithGpsInfoMatrix(:,2), '.b');
% Highlight current track.
hCurTrack = plot(contiPathLossesWithGpsInfo{idxTrack}(:,3), ...
    contiPathLossesWithGpsInfo{idxTrack}(:,2), '.y');
% Highlight the start points of the tracks.
hStartingPoints = plot(contiPathLossesWithGpsInfoFirstPts(:,3), ...
    contiPathLossesWithGpsInfoFirstPts(:,2), 'sg');
% Plot the GPS samples relavant to the signal segment.
hSigSegLocs = plot( ...
    contiPathLossesWithGpsInfo{idxTrack}(indicesRelavantGpsSamps,3), ...
    contiPathLossesWithGpsInfo{idxTrack}(indicesRelavantGpsSamps,2), ...
    '*r');
if exist('axisToSet', 'var')
    axis(axisToSet);
end
plot_google_map('MapType', 'satellite');
legend([hSigSegLocs, hCurTrack, hStartingPoints, hAllGps], ...
    'Signal segment location', 'Current track', ...
    'Starting points', 'All GPS samples', ...
    'Location', 'southeast');
% transparentizeCurLegends;

if exist('fullPathToSavePlot', 'var')
    saveas(hMapFig, fullPathToSavePlot);
    close(hMapFig);
end

end

% EOF