function updateTreeMarkerState(src, evnt)
%UPDATETREEMARKERSTATE Add a new tree marker to markLocs and update the
%overview figure accordingly.
%
% For simplisity, we will just modify variables in the base workspace.
%
% Inputs:
%   - src, evnt
%     The source and event which triggers this callback function.
%
% Yaguang Zhang, Purdue, 05/29/2018

markLocs = evalin('base', 'markLocs');
hInterTreeMarker = evalin('base', 'hInterTreeMarker');
hTreeLocs = evalin('base', 'hTreeLocs');

% (lon, lat).
ptOnClick = evnt.IntersectionPoint(1:2);
% Save (lat, lon).
markLocs = [markLocs; ptOnClick(2) ptOnClick(1)];
disp(['updateTreeMarkerState: New markLoc added!']);
disp(['                       (Lat, Lon) = (', num2str(ptOnClick(2)), ...
    ', ', num2str(ptOnClick(1)), ')']);

% Update the markLocs on the figure.
deleteHandles(hTreeLocs);
[numTrees, ~] = size(markLocs);
hInterTreeMarkerAxes = findall(hInterTreeMarker, 'type', 'axes');
hTreeLocs = plot3(hInterTreeMarkerAxes, ...
    markLocs(:, 2), markLocs(:, 1), ones(numTrees,1), ...
    'rx', 'LineWidth', 2);

% Save the results to the base workspace.
assignin('base', 'markLocs', markLocs);
assignin('base', 'hTreeLocs', hTreeLocs);

end

% EOF