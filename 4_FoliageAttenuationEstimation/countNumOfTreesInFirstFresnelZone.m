function [ numOfTrees ] ...
    = countNumOfTreesInFirstFresnelZone(tx3D, rx3D, trees3D, ...
    fCarrierGHz, FLAG_DEBUG)
%COUNTNUMOFTREESINFIRSTFRESNELZONE Count the number of trees which show up
%in the 1st Fresnel zone.
%
% Input:
%   - tx3D, rx3D
%     Three-element vectors specifying the locations [x y altitude] of the
%     transmitter and reciever, respectively.
%   - trees3D
%     A matrix with each row representing the location of a tree in the
%     form of [x y altitude].
%   - fCarrierGHz
%     The signal frequency in GHz.
%   - FLAG_DEBUG
%     Optional. Set this to be true to show a 2D plot showing all the
%     trees.
% Output:
%   - numOfTrees
%     The number of trees which show up in the 1st Fresnel zone.
%
% For the Fresnel zone radius, we will carry out the computation in 3D;
% Then for determining whether a tree is inside the zone or not, we only
% consider the x-y plane with the radius gotten. Essentially we get the
% distance between the tree and the LOS line in the x-y plane because the
% tree height is not available.
%
% Yaguang Zhang, Purdue, 04/18/2018

if ~exist('FLAG_DEBUG','var')
    FLAG_DEBUG = false;
end

% First, calculate the Fresnel zone radii for the trees.
lambda = physconst( 'LightSpeed' )/(fCarrierGHz*10^9);
trees3DColCell = num2cell(trees3D, 2);
rsFresnel = cellfun(@(p) nFresnelZoneR3D(1, lambda, tx3D, rx3D, p), ...
    trees3DColCell);

% Second, determine whether the trees are within the radius ranges.
losInXY = createLine(tx3D(1:2), rx3D(1:2));
distsTreeToLosInXY = cellfun(@(p) distancePointLine(p(1:2), losInXY), ...
    trees3DColCell);
numOfTrees = sum(distsTreeToLosInXY<rsFresnel);

% Plot for debugging, e.g.
%    [xsT, ysT, zsT] = meshgrid(-200:5:200, -40:5:160, -35:5:105);
%   countNumOfTreesInFirstFresnelZone([-100,100,0], [100,20,70], ...
%       [xsT(:), ysT(:), zsT(:)], 28, true);
if FLAG_DEBUG
    boolsInFirstFresnelZone = distsTreeToLosInXY<rsFresnel;
    ALPHA = 0.2;
    figure; hold on;
    hTx = plot3(tx3D(1), tx3D(2), tx3D(3), '^k');
    hRx = plot3(rx3D(1), rx3D(2), rx3D(3), 'xk');
    hTrees = scatter3(trees3D(~boolsInFirstFresnelZone,1), ...
        trees3D(~boolsInFirstFresnelZone,2), ...
        trees3D(~boolsInFirstFresnelZone,3), '.g', ...
        'MarkerFaceAlpha', ALPHA,'MarkerEdgeAlpha', ALPHA);
    hInZone = scatter3(trees3D(boolsInFirstFresnelZone,1), ...
        trees3D(boolsInFirstFresnelZone,2), ...
        trees3D(boolsInFirstFresnelZone,3), '.b', ...
        'MarkerFaceAlpha', ALPHA,'MarkerEdgeAlpha', ALPHA);
    xlabel('x'); ylabel('y'); zlabel('z');
end

end
% EOF