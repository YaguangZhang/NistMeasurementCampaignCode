function [ catTxt ] = loadCategoryTxt( absPathCatTxt )
%PARSECATEGORYTXT Parse a category .txt file specified by the input absolute
%path.
%
% A category .txt file will have a file name which specifies the category
% (e.g. LoS.txt), and inside it, there is a list of all the measurement
% series belonging to that category (with possibly comment lines indicated
% by a leading '%').
%
% Yaguang Zhang, Purdue, 10/06/2017

% Get the category name.
[~, category] = fileparts(absPathCatTxt);

% Read the string series info in the file.
series = {};
fid = fopen(absPathCatTxt, 'r');
% Read the .txt file line by line.
tline = fgetl(fid);
while ischar(tline)
    if (~strcmp(tline(1),'%'))
        series{end+1} = strtrim(tline); %#ok<AGROW>
    end
    tline = fgetl(fid);
end
fclose(fid);

catTxt = struct();
catTxt.category = category;
catTxt.series = series';

end
% EOF