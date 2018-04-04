function [ ] = deleteHandlesInCell( hs )
%DELETEHANDLESINCELL Check the valibility of handles in the cell hs, and if
%the handle is valid, then delete it.
%
% Yaguang Zhang, Purdue, 11/07/2017

for idx = 1:length(hs)
    if isgraphics(hs{idx}) && isvalid(hs{idx})
        delete(hs{idx});
    end
end

end
% EOF