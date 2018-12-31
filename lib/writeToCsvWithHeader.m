function [ flagSuccess] ...
    = writeToCsvWithHeader(filePath, headerCell, data, precision)
%WRITETOCSVWITHHEADER Write data to a csv file with header.
%
% Inputs:
%   - filePath
%     The full path to the csv file.
%   - headerCell
%     A vector cell representing the elements of the header line.
%   - data
%     A cell or matrix containing the data to be written. 
%   - precision
%     Optional. The precision used for num2str. By default it is '%.8f'.
%
% Yaguang Zhang, Purdue, 12/27/2018

flagSuccess = false;

if ~exist('precision', 'var')
    precision = '%.8f';
end

if iscell(data)
    dataCell = data;
elseif ismatrix(data)
    dataCell = num2cell(data);
end

dataCell = [headerCell; dataCell];

% Write data to end of file
cell2csv(filePath, dataCell, precision);

flagSuccess = true;
end

% EOF