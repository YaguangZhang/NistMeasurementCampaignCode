function [ flagIsS ] = isS( gpggaStr )
%ISW Test whether the location specified by gpggaStr is S (south).
%
% Yaguang Zhang, Purdue, 06/18/2017

strs = strsplit(gpggaStr,',');
flagIsS = strcmp(strs(4), 'S');

end
% EOF