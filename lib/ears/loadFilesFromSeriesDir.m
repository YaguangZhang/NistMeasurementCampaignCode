function [ outFilesDirs, gpsFilesDirs ] ...
    = loadFilesFromSeriesDir( seriesDir, regexPattern )
%LOADFILESFROMSERIESDIR Load the Gnu Radio out files, as well as the
% GPS log files, stored in one _Series folder specified by the input dir
% struct seriesDir.
%
% The outputs, outFilesDirs and gpsFilesDirs, are the resulted dir struct
% arrays, for Gnu Radio out files and the GPS log files, respectively.
%
% By default, all .out files will be loaded. However, if regexPattern for
% the filename is properly set, only the files which agree with the pattern
% will be kept.
%
% Yaguang Zhang, Purdue, 09/26/2017

% Check the folder's name.
[~, seriesFolderName] = fileparts(seriesDir.name);
idxSeriesTokens = regexp(seriesFolderName, ...
    '^Series_(\d+)$', 'tokens');
assert(length(idxSeriesTokens)==1, ...
    'The name for the input folder does not start with `Series_`!');

outFilesDirs = rdir(fullfile(seriesDir.folder, seriesFolderName, ...
    '*.out'), 'regexp(name, ''_\d+.out$'')', true);
if (nargin>1 && ~isempty(regexPattern))
    outFilesDirs = outFilesDirs(arrayfun(@(d) ...
        ~isempty(regexp(d.folder, regexPattern, 'once')), outFilesDirs));
end

gpsFilesDirs = rdir(fullfile(seriesDir.folder, seriesFolderName, ...
    '*_GPS.log'));

end
% EOF