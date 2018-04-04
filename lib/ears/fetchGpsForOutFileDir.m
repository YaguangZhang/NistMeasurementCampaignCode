function [ lat, lon, alt, gpsLog ] ...
    = fetchGpsForOutFileDir( outFileDir )
%FETCHGPSFOROUTFILEDIR Find the cooresponding GPS log file and fetch the
%coordinates stored in it for the .out file specified by the input dir
%struct outFileDir.
%
% Outputs (lat, lon, alt) is the resulted coordinate and gpsLog is the GPS
% information struct containing all the fields stored in the GPS log file.
%
% Yaguang Zhang, Purdue, 09/26/2017

% Get the filename for the .out file.
[~, outFileName] = fileparts(outFileDir.name);

% Find the GPS log file with the same time stamp (specified in the
% filename).
gpsFileDir = rdir(fullfile(outFileDir.folder, [outFileName, '_GPS.log']), '', false);

% Find the closest GPS log as an alternative if on exact match was there.
if isempty(gpsFileDir)
    warning(['No exact GPS log file matched with the timestamp of the input .out file: ', ...
        fullfile(outFileDir.folder, outFileDir.name), ...
        '. We will try to find the GPS log file in this series with the closest timestamp.'])
    
    % Time stamp for the current .out file.
    timeStampOutFile = regexp(outFileName, ...
        'measureSignal_(\d+)$', 'tokens');
    timeStamp = str2num(timeStampOutFile{1}{1});
    
    gpsFilesDirs = rdir(fullfile(outFileDir.folder, '*_GPS.log'), '', false);
    % Get the time stamps.
    timeStampsGpsFiles = arrayfun(@(d) regexp(d.name, ...
    'measureSignal_(\d+)_GPS.log$', 'tokens'), gpsFilesDirs);
    % Convert the stamps (cells) to numbers.
    timeStamps = arrayfun(@(t) str2num(t{1}{1}), timeStampsGpsFiles);
    
    % Find the nearest timestamp.
    [timeDiff, idxGpsFileDir] = min(abs(timeStamps-timeStamp));
    warning(['    GPS file found; Time diff = ', num2str(timeDiff), 's.']);
    
    gpsFileDir = gpsFilesDirs(idxGpsFileDir);
end

% Parse the GPS log.
gpsLog = parseGpsLog(gpsFileDir.name);
gpsLogSample = nmealineread(gpsLog.gpsLocation);

lat = gpsLogSample.latitude;
lon = gpsLogSample.longitude;
alt = gpsLogSample.altitude;

% Add a minus sign if it is W or S.
if(isW(gpsLog.gpsLocation))
    lon = -lon;
end
if(isS(gpsLog.gpsLocation))
    lat = -lat;
end
end
% EOF