function LLR = logLikDifAll(signalMeasurements,timesMean1,timesMean2,Amps1,Amps2,timesNow,sig2,sig_y2)

% Computes the difference in log-likelihoods for two sets of peaks

%%% compute quantities for peak set 1

if isempty(timesMean1)
    dif1 = -signalMeasurements;
else
    % get the signal components
    F1 = signalComponents(timesMean1,timesNow,sig2);
    % get expected signal
    eSig1 = Amps1 * F1';
    % get difference
    dif1 = eSig1 - signalMeasurements;
end

%%% compute quantities for peak set 1

if isempty(timesMean2)
    dif2 = -signalMeasurements;
else
    % get the signal components
    F2 = signalComponents(timesMean2,timesNow,sig2);
    % get expected signal
    eSig2 = Amps2 * F2';
    % get difference
    dif2 = eSig2 - signalMeasurements;
end

%%% get target result

const = -(2*sig_y2)^-1;
llrTermwise = const * (dif1 - dif2) .* (dif1 + dif2);
LLR = sum(llrTermwise,'all');

end