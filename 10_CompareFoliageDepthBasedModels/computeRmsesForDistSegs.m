function [ ituModelPerfTable ] = computeRmsesForDistSegs( ...
    cellDistSegs, allDists, allMeas, cellAllPredicts)
%COMPUTERMSESFORDISTSEGS Compute the RMSE for each distance segment for the
%model predictions provided.
%
%   Inputs:
%       - cellDistSegs
%         A cell with distance segments specified by the min and max values
%         in the form of [min, max), i.e. including min but excluding max.
%       - allDists
%         The distance values we can use for grouping the measurements and
%         predictions accordingly.
%       - allMeas
%         Measurement results.
%       - cellAllPredicts
%         Predictions from different models, each row will be in the form
%         of {'modelName', predictions}.
%   Ouputs:
%       - ituModelPerfTable
%         The results structured as a table.
%
% Yaguang Zhang, Purdue, 09/21/2018

rmseFormatter = '%.2f';

numDistSegs = length(cellDistSegs);
[numModels, ~] = size(cellAllPredicts);

fctRmse = @(measures, predicts) sqrt(mean((predicts - measures).^2));

ituModelPerfCell = cell(numModels, numDistSegs+1);
for idxMod = 1:numModels
    ituModelPerfCell{idxMod, 1} = cellAllPredicts{idxMod, 1};
    for idxDistSeg = 1:numDistSegs
        curDistSegMin = cellDistSegs{idxDistSeg}(1);
        curDistSegMax = cellDistSegs{idxDistSeg}(2);
        
        boolsInCurDistSeg ...
            = (allDists>=curDistSegMin) & (allDists<curDistSegMax);
        
        curPredictions = cellAllPredicts{idxMod, 2}(boolsInCurDistSeg);
        
        ituModelPerfCell{idxMod, 1+idxDistSeg} ...
            = num2str(...
            fctRmse(allMeas(boolsInCurDistSeg), curPredictions), ...
            rmseFormatter);
    end
end

tableHeader = cell(1, numDistSegs+1);
tableHeader{1} = 'ModelName';
for idxDistSeg = 1:numDistSegs
    curDistSegMin = cellDistSegs{idxDistSeg}(1);
    curDistSegMax = cellDistSegs{idxDistSeg}(2);
    
    tableHeader{1+idxDistSeg} = ['Range_', num2str(curDistSegMin), ...
        '_leq_d_lt_', num2str(curDistSegMax)];
end

ituModelPerfTable = cell2table(ituModelPerfCell, ...
    'VariableNames', tableHeader);

end

% EOF