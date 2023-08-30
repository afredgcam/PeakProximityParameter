function log_q  = log_qt_with_prior(signalMeasurements,timesMean,timesNow,amplitudes,peak,eps,sig2)

% This computes the 't' proposal part of the acceptance probability INCLUDING the prior term.

% get indMax
[~,indMax] = size(timesNow);
[~,N] = size(timesMean);

% get the time interval mins and maxs
retainedMeans = [timesMean(1:peak-1),timesMean(peak+1:N)];
zIntMins = max(retainedMeans-eps,1);
zIntMaxs = min(retainedMeans+eps,indMax);

% get which times are out of the zero-regions
regionBool = prod(1 - (timesNow>=zIntMins') .* (timesNow<zIntMaxs'), 1);

% get the times for all q interval ends
[tEnds,tInds] = sort([timesNow,zIntMins,zIntMaxs]); 

% get the (unnormalised) discretised proposal distribution points...
q_ints = timeProbabilities(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2);
q_ints = q_ints .* regionBool;
q = [q_ints,zeros([1,2*(N-1)])];
q = q(tInds);

% ...and normalise it
[~,N_all] = size(q);
avgHeights = 0.5*(q(1:N_all-1) + q(2:N_all));
intervalWidths = tEnds(2:N_all) - tEnds(1:N_all-1);
Q = sum(avgHeights .* intervalWidths);
q = q/Q;

% get probability density of proposing this delay
ti = timesMean(peak);
K_all = find((tEnds > ti) == 1);
K = K_all(1) - 1;
w = (ti - tEnds(K)) / (tEnds(K+1) - tEnds(K));
q_t_new = (1 - w)*q(K) + w*q(K+1);

% get prior ratio term
qNonZero = (q > 0);
intervalsNeeded = qNonZero(2:N_all) | qNonZero(1:N_all-1);
remainingLength = sum(intervalWidths(intervalsNeeded));

% log(q_t_new / 1/rL) = log(q_t_new) + log(rL)
log_q = log(q_t_new + 1e-100) + log(remainingLength);

end