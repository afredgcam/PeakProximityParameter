function log_q  = log_q_anyT(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2)

% This computes the 't' proposal part of the acceptance probability.

% get the (unnormalised) discretised proposal distribution points...
q = timeProbabilities(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2);

% ...and normalise it
N_all = length(q);
avgHeights = 0.5*(q(1:N_all-1) + q(2:N_all));
intervalWidths = timesNow(2:N_all) - timesNow(1:N_all-1);
Q = sum(avgHeights .* intervalWidths);
q = q/Q;

% get probability density of proposing this delay
ti = timesMean(peak);
K_all = find((timesNow > ti) == 1);
K = K_all(1) - 1;
w = (ti - timesNow(K)) / (timesNow(K+1) - timesNow(K));
q_t_new = (1 - w)*q(K) + w*q(K+1);

% log(q_t_new)
log_q = log(q_t_new + 1e-100);

end