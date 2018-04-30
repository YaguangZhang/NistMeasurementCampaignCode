function [ allSigOutFiles ] = locateFilteredOutFiles(absPathToSearch)
%LOCATEFILTEREDOUTFILES Search for all the filtered .out
%signal recording files recursively under absPathToSearch.
%
% The output allSigOutFiles is a structure array created by rdir.
%
% Yaguang Zhang, Purdue, 04/30/2018

% Find all GnuRadio .out files.
allSigOutFiles = rdir(fullfile(absPathToSearch, '**', '*.out'));
% We do not need the filtered version.
allSigOutFiles = allSigOutFiles(arrayfun(@(p) ...
    ~isempty(regexp(p.name, '\d+_filtered.out','match')), allSigOutFiles));

end
% EOF
