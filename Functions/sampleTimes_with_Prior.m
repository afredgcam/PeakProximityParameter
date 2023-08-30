function [t_i, log_q]  = sampleTimes_with_Prior(signalMeasurements,timesMean,timesNow,amplitudes,peak,eps,sig2,wanting_q)

% This samples the mean of one peak and computes part of the acceptance probability INCLUDING the prior term.

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

% ...and the CDF of the continuous version at those same points
[~,N_all] = size(q);
avgHeights = 0.5*(q(1:N_all-1) + q(2:N_all));
intervalWidths = tEnds(2:N_all) - tEnds(1:N_all-1);
Q = zeros([1,N_all]);
Q(2:N_all) = cumsum(avgHeights .* intervalWidths);

% normalise them
q = q/Q(N_all);
Q = Q/Q(N_all);

% get uniform sample
u = rand(1);

% find the section of continuous proposal 
K_all = find((Q > u) == 1);
K = K_all(1) - 1;

% get sampled ti
q_diff = q(K+1) - q(K);
T_diff = tEnds(K+1) - tEnds(K);
if q_diff ~= 0
    radical_squared = q(K)^2 + 2 * q_diff/T_diff * (u - Q(K));
    ti = tEnds(K) - T_diff/q_diff * (q(K) - radical_squared^0.5);
elseif q(K) ~= 0
    ti = tEnds(K) + (u - Q(K))/q(K);
end

% ensure that it is in the range, and not at the very end points
ti = min([ti,indMax-0.0001]);
ti = max([ti,1.0001]);

% assemble the current times
t_i = timesMean;
t_i(peak) = ti;

if wanting_q
    % get probability density of proposing this ti
    w = (ti - tEnds(K)) / (tEnds(K+1) - tEnds(K));
    q_t_new = (1 - w)*q(K) + w*q(K+1);

    % get prior ratio term
    qNonZero = (q > 0);
    intervalsNeeded = qNonZero(2:N_all) | qNonZero(1:N_all-1);
    remainingLength = sum(intervalWidths(intervalsNeeded));

    % log(q_t_new / 1/rL)
    log_q = log(q_t_new + 1e-100) + log(remainingLength);
else
    log_q = 0;
end

end