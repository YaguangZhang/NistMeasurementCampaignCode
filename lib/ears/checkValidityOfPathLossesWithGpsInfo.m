function [ boolsValidPathlosses, ...
    boolsInvalidCoor, boolsInfPathloss, boolsInvalidData ] ...
    = checkValidityOfPathLossesWithGpsInfo(pathLossesWithGpsInfo, ...
    relPathsOutFilesUnderDataFolder)
%CHECKVALIDITYOFPATHLOSSESWITHGPSINFO Check the validity of the oupts
%pathLossesWithGpsInfo from evalPathLosses.
%
% Output:
%    - boolsValidPathlosses
%      A bool vector indicating whether the row records in
%      pathLossesWithGpsInfo are valid.
%    - relPathsOutFilesUnderDataFolder
%      Optional. When present, it contains the relative paths for the row
%      records in boolsValidPathlosses, and we will check them for invalid
%      data listed in relPathSegInvalidData.
%
% Yaguang Zhang, Purdue, 10/16/2017

if nargin >1
    % Some data collected are not valid from the very begining. For example,
    % the '20170619_LargeScale/Series_9' data are not valid because the
    % measurement was interrupted by a rain.
    %   Todo: add supports for multiple directories.
    relPathSegInvalidData = fullfile('20170619_LargeScale', 'Series_9');
    
    boolsInvalidData = cellfun(@(p) contains(p, relPathSegInvalidData), ...
        relPathsOutFilesUnderDataFolder);
    if any(boolsInvalidData)
        warning([num2str(sum(boolsInvalidData)), ...
            ' specified invalid data detected.', ...
            ' We will ignore the cooresponding results.']);
    end
else
    [numRowRecords, ~] = size(pathLossesWithGpsInfo);
    boolsInvalidData = zeros(numRowRecords, 1);
end

% Then, check the GPS info.
disp('checkValidityOfPathLossesWithGpsInfo ...');
boolsInvalidCoor = pathLossesWithGpsInfo(:,2)==0 ...
    & pathLossesWithGpsInfo(:,3)==0;
if any(boolsInvalidCoor)
    warning([num2str(sum(boolsInvalidCoor)), ...
        ' invalid (lat, lon) pairs detected (both are 0).', ...
        ' We will ignore these points together with their path losses.']);
end

% At last, the inf path losses.
boolsInfPathloss = isinf(pathLossesWithGpsInfo(:,1));
if any(boolsInfPathloss)
    warning([num2str(sum(boolsInfPathloss)), ...
        ' inf path loss detected.', ...
        ' We will show these points at z=0 with different markers.']);
end

boolsValidPathlosses = ...
    ~(boolsInvalidCoor | boolsInfPathloss | boolsInvalidData);

disp('Done!');
end
% EOF