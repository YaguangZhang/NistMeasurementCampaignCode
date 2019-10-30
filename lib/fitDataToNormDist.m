function [ADTestResult, ADTestP, muHat, sigmaHat, hFig] ...
    = fitDataToNormDist(x)
%FITDATATONORMDIST Carry out an Aderson-Darling test on the input data and
%fit them to a norm distribution.
%
% Input:
%   - x
%     The input data as a vector.
% Outputs:
%   - ADTestResult
%     The test decision for the null hypothesis (0) that the data in vector
%     x is from a population with a normal distribution, using the
%     Anderson-Darling test. The alternative hypothesis is that x is not
%     from a population with a normal distribution (1). The result h is 1
%     if the test rejects the null hypothesis at the 5% significance level,
%     or 0 otherwise.
%   - ADTestP
%     p-value of the Anderson-Darling test, returned as a scalar value in
%     the range [0,1]. ADTestP is the probability of observing a test
%     statistic as extreme as, or more extreme than, the observed value
%     under the null hypothesis.
%   - muHat, sigmaHat
%     Estimates of normal distribution parameters (the mean muHat and
%     standard deviation sigmaHat).
%   - hFig
%     Optional output for the handler of the illustration figure. The
%     illustration figure for the results will be generated only if this
%     output is present.
%
% Yaguang Zhang, Purdue, 10/16/2019

% By default, find the flag flagGenFigSilently in the base workspace. If
% not found, flagGenFigSilently will be set to false.
try
    flagGenFigSilently = evalin('base', 'flagGenFigSilently');
catch
    flagGenFigSilently = false;
end

MIN_NUM_OF_BINS = 20;

if isempty(x)
    [ADTestResult, ADTestP, muHat, sigmaHat] = deal(nan);
    hFig = figure('visible', ~flagGenFigSilently);
else
    [ADTestResult, ADTestP] = adtest(x);
    [muHat, sigmaHat] = normfit(x);
    
    if nargout>4
        % For plotting the impirical pdf.
        numBins = max(ceil(length(x)./10), MIN_NUM_OF_BINS);
        [N, edges] = histcounts(x, numBins);
        empiricalPdfXs = mean([edges(1:end-1); edges(2:end)]);
        empiricalPdfYs = N./length(x)./(edges(2)-edges(1));
        
        % For plotting the fitted norm dist.
        numNormPts = 100;
        xsNorm = linspace(edges(1), edges(end), numNormPts);
        ysNorm = normpdf(xsNorm, muHat, sigmaHat);
        
        % For showing extra statistics in the plot title.
        skew = skewness(x);
        skewUnbiased = skewness(x, 0);
        kurt = kurtosis(x);
        kurtUnbiased = kurtosis(x, 0);
        mu = mean(x);
        sigma = std(x);
        
        % Float number fommatters.
        FNF1Pre = '%.1f';
        FNF2Pre = '%.2f';
        
        % Figure.
        hFig = figure('visible', ~flagGenFigSilently); hold on;
        hEmp = plot(empiricalPdfXs, empiricalPdfYs, '--b');
        hNorm = plot(xsNorm, ysNorm, '--r');
        axis tight; grid on; grid minor;
        xlabel('Value'); ylabel('PDF');
        legend([hEmp, hNorm], 'Empirical', 'Fitted norm')
        title({['mu = ', num2str(mu, FNF1Pre), ...
            ', sigma = ', num2str(sigma, FNF1Pre), ...
            ', Anderson-Darling Test Result: ', num2str(ADTestResult), ...
            ' (p = ', num2str(ADTestP, FNF2Pre), ')']; ...
            ['Skewness = ', num2str(skew, FNF2Pre), ...
            ' (unbiased:', num2str(skewUnbiased, FNF2Pre), '), ', ...
            'Kurtosis = ', num2str(kurt, FNF2Pre), ...
            ' (unbiased:', num2str(kurtUnbiased, FNF2Pre), ')']; ...
            ['Fitted Norm: muHat = ', num2str(muHat, FNF1Pre), ...
            ', sigmaHat = ', num2str(sigmaHat, FNF1Pre)]});
    end
end
end
% EOF