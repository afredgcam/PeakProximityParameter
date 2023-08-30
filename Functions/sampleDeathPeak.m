function [peak,log_q_peak] = sampleDeathPeak(A)

% get peak probabilities
log_probs = peakDeathProb(A);
probs = exp(log_probs);

% sample peak
N = length(A);
peak = randsample(N,1,true,probs);

% save log(prob) of that peak's selection
log_q_peak = log_probs(peak);

end