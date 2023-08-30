function F = signalComponents(timesMean,timesNow,variance)

% Computes the F matrix, used for computing the expected echo curve

% Inputs:
    % timesMean ((1 x N) array): some sampled means, {µ}
    % timesNow ((1 x K) array): the sampling times, t
    % variance (float): the peak variance, σ2

% Outputs:
    % F ((K x N) array): the desired F matrix



[~,numSent] = size(timesMean); % number of pulses (i.e., N)
[~,numNow] = size(timesNow); % number of indices

% make sizes compatible for 'normpdf()'
Sent = repmat(timesMean,numNow,1);
Now = repmat(timesNow',1,numSent);

% get desired matrix result
F = normpdf(Now,Sent,variance^0.5); % it's a noNow x N matrix

end