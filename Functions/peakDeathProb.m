function log_probs = peakDeathProb(A)

% get current number of peaks
N = length(A);

% get log(probabities)
log_probs = log(abs(A).^-1) - log(sum(abs(A).^-1));

end