function [latM, lonM, altM] = fetchMedianGpsForSeriesDir( absPathSeriesDir )
%FETCHMEDIANGPSFORSERIESDIR Find all the cooresponding GPS log files, fetch
%the coordinates stored in them, and average the locked GPS samples via
%median, for the Series folder specified by the input absolute path.
%
% Yaguang Zhang, Purdue, 10/04/2017

% Find all the GPS log files.
gpsLogsDirs = rdir(fullfile(absPathSeriesDir, '**', '*_GPS.log'), ...
    '', false);
% Keep a record of (lat, lon, alt, locked).
gpsRecords = nan(length(gpsLogsDirs), 4);

for idxGpsLog = 1:length(gpsLogsDirs)
    % Parse the GPS log.
    gpsLog = ...
        parseGpsLog(gpsLogsDirs(idxGpsLog).name);
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
    
    locked = str2double(gpsLog.gpsLocked);
    
    % Store the result.
    gpsRecords(idxGpsLog,:) = [lat, lon, alt, locked]; 
end

gpsRecordsLocked = gpsRecords(gpsRecords(:,4)==1,:);

latM = median(gpsRecordsLocked(:,1));
lonM = median(gpsRecordsLocked(:,2));
altM = median(gpsRecordsLocked(:,3));

end

% EOF