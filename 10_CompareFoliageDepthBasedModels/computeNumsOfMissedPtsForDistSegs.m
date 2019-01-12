function [ modelPerfTableNumsOfMissedPts, modelPerfCellNumsOfMissedPts, ...
    modelPerfTableAccuracies, modelPerfCellAccuracies] ...
    = computeNumsOfMissedPtsForDistSegs( ...
    cellDistSegs, allDists, allMeas, cellAllPredicts, maxAbsErrAllowedInDb)
%COMPUTE SFORDISTSEGS Compute the number of missed prediction points (and
%the according accuracy values) for each distance segment for the model
%predictions provided.
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
%       - maxAbsErrAllowedInDb
%         The threshold for considering a predicaiton value to be valid or
%         not. Values with absolute errors over this parameter will be
%         considered as missed prediction points.
%   Ouputs:
%       - modelPerfTable
%         The results structured as a table.
%       - modelPerfCell
%         The results structured as a table.
%
% Yaguang Zhang, Purdue, 09/21/2018

numsOfMissedPtsFormatter = '%d';
accuraciesFormatter = '%.2f';

numDistSegs = length(cellDistSegs);
[numModels, ~] = size(cellAllPredicts);

fctNumsOfMissedPts = @(measures, predicts) ...
    sum(abs(predicts - measures)>maxAbsErrAllowedInDb);

modelPerfCellNumsOfMissedPts = cell(numModels, numDistSegs+1);
modelPerfStrCellNumsOfMissedPts = cell(numModels, numDistSegs+1);
for idxMod = 1:numModels
    modelPerfCellNumsOfMissedPts{idxMod, 1} = cellAllPredicts{idxMod, 1};
    modelPerfStrCellNumsOfMissedPts{idxMod, 1} = cellAllPredicts{idxMod, 1};
    
    for idxDistSeg = 1:numDistSegs
        curDistSegMin = cellDistSegs{idxDistSeg}(1);
        curDistSegMax = cellDistSegs{idxDistSeg}(2);
        
        boolsInCurDistSeg ...
            = (allDists>=curDistSegMin) & (allDists<curDistSegMax);
        
        curPredictions = cellAllPredicts{idxMod, 2}(boolsInCurDistSeg);        
        curNumOfMissedPts ...
            = fctNumsOfMissedPts(allMeas(boolsInCurDistSeg), ...
            curPredictions);
        
        modelPerfCellNumsOfMissedPts{idxMod, 1+idxDistSeg} = curNumOfMissedPts;
        modelPerfStrCellNumsOfMissedPts{idxMod, 1+idxDistSeg} ...
            = num2str(curNumOfMissedPts, numsOfMissedPtsFormatter);
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

modelPerfTableNumsOfMissedPts = cell2table(modelPerfStrCellNumsOfMissedPts, ...
    'VariableNames', tableHeader);

modelPerfCellAccuracies = cell(numModels, numDistSegs+1);
modelPerfStrCellAccuracies = cell(numModels, numDistSegs+1);
for idxMod = 1:numModels
    modelPerfCellAccuracies{idxMod, 1} = cellAllPredicts{idxMod, 1};
    modelPerfStrCellAccuracies{idxMod, 1} = cellAllPredicts{idxMod, 1};
    
    for idxDistSeg = 1:numDistSegs
        curDistSegMin = cellDistSegs{idxDistSeg}(1);
        curDistSegMax = cellDistSegs{idxDistSeg}(2);
        
        boolsInCurDistSeg ...
            = (allDists>=curDistSegMin) & (allDists<curDistSegMax);
        curTotalNumPts = sum(boolsInCurDistSeg);
        
        curPredictions = cellAllPredicts{idxMod, 2}(boolsInCurDistSeg);
        curAccuracy = 1-fctNumsOfMissedPts(allMeas(boolsInCurDistSeg), ...
            curPredictions)/curTotalNumPts;
        
        modelPerfCellAccuracies{idxMod, 1+idxDistSeg} = curAccuracy;
        modelPerfStrCellAccuracies{idxMod, 1+idxDistSeg} ...
            = num2str(curAccuracy, accuraciesFormatter);
    end
end

modelPerfTableAccuracies = cell2table(modelPerfStrCellAccuracies, ...
    'VariableNames', tableHeader);

end

% EOF