function x_hat = averageCurve(meanSamples,amplitudeSamples,interval,sig2)

% get number of data points / number of samples
K = length(interval);
S = length(meanSamples);

% make all curves
x_all = zeros([S,K]);
for sample = 1:S
    if ~isempty(meanSamples{sample})
        x_all(sample,:) = amplitudeSamples{sample} * signalComponents(meanSamples{sample},interval,sig2)';
    end
end

% get the avereage over all curves
x_hat = mean(x_all,1);

end