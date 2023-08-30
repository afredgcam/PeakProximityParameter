function probT = timeProbabilities(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2)

% Computes proposal density at sampling times (that is, q( tk | … ) )

% Inputs:
    % signalMeasurements ((1 x K) array): the echo curve, y
    % timesMean ((1 x N) array): the most recent mean samples, {µ1,...,µ{i-1}}j
        % ⋃ {µi,...,µN}{j-1}
    % timesNow ((1 x K) array): the sampling times, t
    % amplitudes ((1 x N) array): the current amplitude samples, Aj
    % peak (integer): the peak that needs resampling, i
    % sig2 (float): the peak variance, σ2

% Outputs:
    % probT ((1 x K) array): proposal density at sampling times, q(tk|…)


% get number of peaks
[~,N] = size(timesMean);

% remove the current peak
tProxy = [timesMean(1:peak-1),timesMean(peak+1:N)];
AProxy = [amplitudes(1:peak-1),amplitudes(peak+1:N)];

% get expected signal from these components
expSignal = signalComponents(tProxy,timesNow,sig2) * AProxy';

% get unexplained data with these peaks
unexData = signalMeasurements - expSignal';

% discretise gaussian peak with same amp as original (only sign matters) <-- not sure sign matters anymore because I now take absolute values...
singlePeak = normpdf(-3:sig2^-0.5:3) * (2*(amplitudes(peak) > 0) - 1);

% get convolution
probT = abs(conv(unexData,singlePeak,'same'));

end