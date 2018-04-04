function [ gpsSample ] = parseGpsLog( filePath )
%PARSEGPSLOG Summary of this function goes here
%   Parse the GPS log file generated in the EARS measurement campaign.
%
% Input: 
%   - filePath
%     The full path to the GPS log file. Supports csv type of string rows
%     in the file (Note: ", " with 1 space). All "{" and "}" will be
%     removed.
% Output: 
%   - gpsSample
%     A struct with fields generated according to the title line of the
%     file. Each field contains the corresponding value as a string. 
%
% Yaguang Zhang, Purdue, 06/18/2017

% Read in the first line as fields and the second one as values.
fid = fopen(filePath);
fields = strsplit(regexprep(fgetl(fid),'[{}]',''), ', ');
values = strsplit(regexprep(fgetl(fid),'[{}]',''), ', ');
fclose(fid);

% Output struct.
gpsSample = cell2struct(values,fields,2);

end
% EOF