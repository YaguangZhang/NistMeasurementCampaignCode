function [lat, lon, alt, gpsLogSample] = parseNmeaStr(gpsLocation)
%PARSENMEASTR Parse the NMEA GPS string.
%
% Yaguang Zhang, Purdue, 04/09/2018

if isempty(gpsLocation)
    % The NMEA string is empty.
    [lat, lon, alt, gpsLogSample] = deal(nan);
else    
    gpsLogSample = nmealineread(gpsLocation);
    
    lat = gpsLogSample.latitude;
    lon = gpsLogSample.longitude;
    alt = gpsLogSample.altitude;
    
    % Add a minus sign if it is W or S.
    if(isW(gpsLocation))
        lon = -lon;
    end
    if(isS(gpsLocation))
        lat = -lat;
    end
end

end
% EOF