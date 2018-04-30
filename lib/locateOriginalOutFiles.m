function [ allSigOutFiles ] = locateOriginalOutFiles(absPathToSearch)
%LOCATEORIGINALOUTFILES Search for all the original (non-filtered) .out
%signal recording files recursively under absPathToSearch.
%
% The output allSigOutFiles is a structure array created by rdir.
%
% Yaguang Zhang, Purdue, 04/30/2018

% Find all GnuRadio .out files.
allSigOutFiles = rdir(fullfile(absPathToSearch, '**', '*.out'));
% We do not need the filtered version.
allSigOutFiles = allSigOutFiles(arrayfun(@(p) ...
    ~isempty(regexp(p.name, '\d+.out','match')), allSigOutFiles));

end
% EOF