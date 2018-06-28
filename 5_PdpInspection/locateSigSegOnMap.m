function [ hMapFig ] = locateSigSegOnMap(...
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
%   Output:
%       - hMapFig
%         The map plot generated to show the singal segment's location.
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

hMapFig = figure; hold on;
% Plot all the GPS samples.
hAllGps = plot(contiPathLossesWithGpsInfoMatrix(:,3), ...
    contiPathLossesWithGpsInfoMatrix(:,2), '.b');
% Plot the GPS samples relavant to the signal segment.
hSigSegLocs = plot( ...
    contiPathLossesWithGpsInfo{idxTrack}(indicesRelavantGpsSamps,3), ...
    contiPathLossesWithGpsInfo{idxTrack}(indicesRelavantGpsSamps,2), ...
    '*r');
if exist('axisToSet', 'var')
    axis(axisToSet);
end
plot_google_map('MapType', 'satellite');
legend([hSigSegLocs, hAllGps], ...
    'Signal segment location', 'All GPS samples', ...
    'Location', 'southeast');
% transparentizeCurLegends;

if exist('fullPathToSavePlot', 'var')
    saveas(hMapFig, fullPathToSavePlot);
    close(hMapFig);
end

end

% EOF