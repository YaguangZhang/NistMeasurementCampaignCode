function saveMarkLocs(~, ~)
%SAVEMARKLOCS Save the variable markLocs.
%
%   We will save the variable markLocs in the base workspace into the file
%   specified by the variable ABS_PATH_TO_SAVE_TREE_LOCS in the base
%   workspace.
%
% Inputs:
%   - src, evnt
%     The source and event which triggers this callback function.
%
% Yaguang Zhang, Purdue, 05/29/2018

markLocs = evalin('base', 'markLocs');
ABS_PATH_TO_SAVE_TREE_LOCS = evalin('base', 'ABS_PATH_TO_SAVE_TREE_LOCS');

save(ABS_PATH_TO_SAVE_TREE_LOCS, 'markLocs');
disp('saveMarkLocs: markLocs saved to file!')

closereq;

end

% EOF