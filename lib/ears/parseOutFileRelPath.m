function [ date, type, serNum, timestamp ] ...
    = parseOutFileRelPath( relPathOutFile )
%PARSEOUTFILERELPATH Parse a relative path for a .out file.
% 
% Input:
%    - relPathOutFile
%      A string to specify the relative path for the .out file, e.g.
%        '20170617_LargeScale/Series_1/measureSignal_1497709963.out'
%      or
%        '20170617_LargeScale\Series_1\measureSignal_1497709963.out'
%
% Outputs:
%    - date
%      Date string, e.g. '20170617'.
%    - type
%      Measurement type string, 'LargeScale, 'SIMO', or 'Conti'.
%    - serNum
%      A scalar for the series number, e.g. 1.
%    - timestamp
%      A scalar for the timestamp, e.g. 1497709963.
%
% Yaguang Zhang, Purdue, 10/16/2017

tokens = regexp(relPathOutFile, ...
    '(\d+)_([a-zA-Z]+)[\\\/]Series_(\d+)[\\\/]measureSignal_(\d+).out$', ...
    'tokens');
tokens = tokens{1};
date = tokens{1, 1};
type = tokens{1, 2};
serNum = str2num(tokens{1, 3});
timestamp = str2num(tokens{1, 4});

end
% EOF