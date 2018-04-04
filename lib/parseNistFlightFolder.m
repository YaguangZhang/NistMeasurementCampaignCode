function [ samples, globalSettings ] = parseNistFlightFolder( absPathToFlightFolder, ...
    secondsFromEnd )
%PARSENISTFLIGHTFOLDER Load the data under one flight folder.
%
% V1: We only care about the USRP signal recording.
%
% Inputs:
%   - absPathToFlightFolder
%     A string specifying the absolute path for the flight folder to parse.
%     Example folder name: 'flight_number_15'.
%   - secondsFromEnd
%     Optional. An integer specifying how many seconds will be loaded from
%     the end of the recording. If not present, all samples will be loaded.
%
% Outputs:
%   - samples
%     The loaded complex signal.
%   - globalSettings
%     The parsed global settings recorded in global.json.
%
% Yaguang Zhang, Purdue, 03/30/2018

% JSON file for global settings.
absPathToGlobal = fullfile(absPathToFlightFolder, 'global.json');
% Binary file for the complex USPR recording.
absPathToSigFile = fullfile(absPathToFlightFolder, 'samples.sigmf-data');

% Read in and parse the global settings.
fGlobal = fopen(absPathToGlobal);
strGlobal = char(fread(fGlobal,inf)');
fclose(fGlobal);
globalSettings = jsondecode(strGlobal);

% Read in the signal.
if nargin>1
    count = secondsFromEnd.*globalSettings.core_sample_rate;
    samples = readComplexBinaryFromEnd(absPathToSigFile, count);
else
    samples = read_complex_binary(absPathToSigFile);
end

end
% EOF