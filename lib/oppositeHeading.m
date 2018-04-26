function [ oppHeading ] = oppositeHeading( heading )
%OPPOSITEHEADING Compute the opposite heading in true-north east manner.
%
% Input:
%   - heading
%     Matrix with elements from (0, 360] in degree.
%
% Output:
%   - oppHeading
%     Matrix with elements from (0, 360] in degree. They coorespond to the
%     opposite directions specified by the elements of heading.
%
% Yaguang Zhang, Purdue, 05/17/2017

oppHeading = mod(heading+180, 360);
oppHeading(oppHeading==0) = 360;

end
% EOF