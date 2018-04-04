function [indicesStarts, indicesEnds, values] ...
    = findConsecutiveSubSeqs(parentSeq)
% FINDCONSECUTIVESUBSEQS Find consecutive subsequences.
%
% Inputs:
%
%   - parentSeq
%
%   The sequence "parentSeq" which will be searched for subsequences of
%   consecutive elements.
%
% Outputs:
%
%   - indicesStarts, indicesEnds, values
%
%   The start indices, end indices and the values for those subsequences
%   found.
%
% Yaguang Zhang, Purdue, 09/15/2015

diffParentSeq = diff(parentSeq);

% Find consecutive 0's in the diff vector.
[indicesStarts, indicesEnds] = findConsecutiveSubSeq(diffParentSeq, 0);
indicesEnds = indicesEnds + 1;
 
values = parentSeq(indicesStarts);

end
% EOF