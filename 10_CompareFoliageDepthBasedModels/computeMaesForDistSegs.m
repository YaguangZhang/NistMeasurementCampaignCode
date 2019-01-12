function [ modelPerfTable, modelPerfCell ] ...
    = computeMaesForDistSegs( ...
    cellDistSegs, allDists, allMeas, cellAllPredicts)
%COMPUTEMAESFORDISTSEGS Compute the Mean Absolute Error (MAE) for each
%distance segment for the model predictions provided.
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
%       - modelPerfTable
%         The results structured as a table.
%       - modelPerfCell
%         The results structured as a table.
%
% Yaguang Zhang, Purdue, 09/21/2018

maeFormatter = '%.2f';

numDistSegs = length(cellDistSegs);
[numModels, ~] = size(cellAllPredicts);

fctMae = @(measures, predicts) mean(abs(predicts - measures));

modelPerfCell = cell(numModels, numDistSegs+1);
modelPerfStrCell = cell(numModels, numDistSegs+1);
for idxMod = 1:numModels
    modelPerfCell{idxMod, 1} = cellAllPredicts{idxMod, 1};
    modelPerfStrCell{idxMod, 1} = cellAllPredicts{idxMod, 1};
    
    for idxDistSeg = 1:numDistSegs
        curDistSegMin = cellDistSegs{idxDistSeg}(1);
        curDistSegMax = cellDistSegs{idxDistSeg}(2);
        
        boolsInCurDistSeg ...
            = (allDists>=curDistSegMin) & (allDists<curDistSegMax);
        
        curPredictions = cellAllPredicts{idxMod, 2}(boolsInCurDistSeg);        
        curMae = fctMae(allMeas(boolsInCurDistSeg), curPredictions);
        
        modelPerfCell{idxMod, 1+idxDistSeg} = curMae;
        modelPerfStrCell{idxMod, 1+idxDistSeg} ...
            = num2str(curMae, maeFormatter);
    end
end

tableHeader = cell(1, numDistSegs+1);
tableHeader{1} = 'ModelName';
for idxDistSeg = 1:numDistSegs
    curDistSegMin = cellDistSegs{idxDistSeg}(1);
    curDistSegMax = cellDistSegs{idxDistSeg}(2);
    
    tableHeader{1+idxDistSeg} ...
        = ['Range_', strrep(num2str(curDistSegMin), '.', '_'), ...
        '_leq_d_lt_', strrep(num2str(curDistSegMax), '.', '_')];
end

modelPerfTable = cell2table(modelPerfStrCell, ...
    'VariableNames', tableHeader);

end

% EOF