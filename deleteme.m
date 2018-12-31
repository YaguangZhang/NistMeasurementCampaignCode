for m = 1:10%numRows
    if (mod(m, floor(numRows/100))==0)
        disp(['    Row ', num2str(m), '/' num2str(numRows)]);
    end
    dataRowCnt = 1;
    for n = 1:numCols
        curData(dataRowCnt,:) = {m, n, ...
            VEG_AREA_IMG_META.XS(m,n), VEG_AREA_IMG_META.YS(m,n), ...
            VEG_AREA_IMG_META.UTM_ZONE, ...
            VEG_AREA_IMG_META.LATS(m,n), VEG_AREA_IMG_META.LONS(m,n), ...
            VEG_AREA_IMG_META.ALTS(m,n), ...
            groundHeightWrtTXInM(m,n), treeHeightInM(m,n), ...
            vegAreas(m,n)};
        dataRowCnt = dataRowCnt+1;
    end
    disp('  Writing a batch to file ...')    
    appendCell2Csv(fullPathFoliageAreaCsv, curData, '%.8f');
    disp('  Done!')
end
disp('Done!')