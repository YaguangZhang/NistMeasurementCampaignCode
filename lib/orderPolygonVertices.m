function [orderedXs, orderedYs, orderedIndices] ...
    = orderPolygonVertices(polygonXs, polygonYs)
%ORDERPOLYGONVERTICES Order the input cortices for a polygon
%counter-clock-wise to simply the shape.
%
% Inputs:
%   - polygonXs, polygonXs
%     Column vectors for the x and y coordinates of the input polygon,
%     respectively.
%
% Outputs:
%   - orderedXs, orderedYs
%     Column vectors for theColumn vectors for the x and y coordinates
%     output ordered x and y coordinates, respectively.
%   - orderedIndices
%     The indices to indicate the new order of the input coordinates to
%     form the output coordinates.
%
% Yaguang Zhang, Purdue, 10/14/2019

assert((polygonXs(1)==polygonXs(end))&&(polygonYs(1)==polygonYs(end)), ...
    'The input polygon should be closed!');

% Remove the last repeated point.
polygonXs(end) = [];
polygonYs(end) = [];

polygonCentroidX = mean(polygonXs);
polygonCentroidY = mean(polygonYs);
[~, orderedIndices] = sort(atan2(polygonXs-polygonCentroidX, ...
    polygonYs-polygonCentroidY));

% Close the output polygon.
orderedIndices(end+1) = orderedIndices(1);

orderedXs = polygonXs(orderedIndices);
orderedYs = polygonYs(orderedIndices);

end
% EOF