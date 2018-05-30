function [] = deleteHandles(handles)
%DELETEHANDLES A helper to delete handles in a cell/vector.
%
% Yaguang Zhang, Purdue, 08/26/2015

if iscell(handles)
    % Assume each element of the cell is a handle/ a handle vector.
    for idxHandle = 1:numel(handles)
        if ~isempty(handles{idxHandle})
            if isvalid(handles{idxHandle})
                delete(handles{idxHandle});
            end
        end
    end
else
    if isvalid(handles)
        % Assume the input is a handle/ a handle vector.
        delete(handles);
    end
end

end

% EOF