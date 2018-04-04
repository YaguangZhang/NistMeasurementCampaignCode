function [indicesStarts, indicesEnds] = findConsecutiveSubSeq(parentSeq, element)
% FINDCONSECUTIVESUBSEQ Find consecutive subsequences.
%
% Inputs:
%
%   - parentSeq, element
%
%   The sequence "parentSeq" which will be searched for subsequences of
%   consecutive "element".
%
% Outputs:
%
%   - indicesStarts, indicesEnds
%
%   Column vectors for the start indices and end indices, respectively, of
%   those subsequences found.
%
% Yaguang Zhang, Purdue, 03/11/2015

% We will use here similar code for finding -100 sequences before. In
% order to do that, we need to change "element" to -1 and all
% other labels (higher densities) to 0.
diffParentSeq = -ones(length(parentSeq),1);
diffParentSeq(parentSeq ~= element) = 0;
diffParentSeq = diff([0;diffParentSeq;0]);

% In terms of 1:indicesUnlabeledSeq.
indicesStarts = find(diffParentSeq==-1); % start
indicesEnds = find(diffParentSeq==1)-1; % end

end
% EOF