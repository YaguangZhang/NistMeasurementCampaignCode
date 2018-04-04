function hOutFig = cloneFig(hInFig,hOutFig)
%CLONEFIG Clone one figure to another.
%
% Modified from CloneFig from:
%   - Matt Fetterman, 2009 pretty much taken from Matlab Technical
%   solutions:
% http://www.mathworks.com/support/solutions/en/data/1-1UTBOL/?solution=1-1UTBOL
%
% Yaguang Zhang, Purdue, 09/13/2017

if nargin<2
    hOutFig = figure;
end
set(0, 'CurrentFigure', hOutFig);
clf;
compCopy(hInFig,hOutFig);

function compCopy(op, np)
%COMPCOPY copies a figure object represented by "op" and its descendants to
%another figure "np" preserving the same hierarchy.
ch = get(op, 'children');
if ~isempty(ch)
    nh = copyobj(ch,np);
    for k = 1:length(ch)
        compCopy(ch(k),nh(k));
    end
end
return
%EOF