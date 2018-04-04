function [ flagIsW ] = isW( gpggaStr )
%ISW Test whether the location specified by gpggaStr is W (west).
%
% Yaguang Zhang, Purdue, 06/18/2017

strs = strsplit(gpggaStr,',');
flagIsW = strcmp(strs(6), 'W');

end
% EOF