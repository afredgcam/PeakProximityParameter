function [t_i, log_q]  = sampleAnyTime(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2,wanting_q)

% This samples the mean of one peak and computes part of the acceptance probability.

% get the (unnormalised) discretised proposal distribution points...
q = timeProbabilities(signalMeasurements,timesMean,timesNow,amplitudes,peak,sig2);

% ...and the CDF of the continuous version at those same points
[~,N_all] = size(q);
avgHeights = 0.5*(q(1:N_all-1) + q(2:N_all));
intervalWidths = timesNow(2:N_all) - timesNow(1:N_all-1);
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
T_diff = timesNow(K+1) - timesNow(K);
if q_diff ~= 0
    radical_squared = q(K)^2 + 2 * q_diff/T_diff * (u - Q(K));
    ti = timesNow(K) - T_diff/q_diff * (q(K) - radical_squared^0.5);
elseif q(K) ~= 0
    ti = timesNow(K) + (u - Q(K))/q(K);
end

% ensure that it is in the range, and not at the very end points
ti = min([ti,timesNow(N_all)-0.01]);
ti = max([ti,timesNow(1)+0.01]);

% assemble the current times
t_i = timesMean;
t_i(peak) = ti;

if wanting_q
    % get probability density of proposing this ti
    w = (ti - timesNow(K)) / (timesNow(K+1) - timesNow(K));
    q_t_new = (1 - w)*q(K) + w*q(K+1);

    % log(q_t_new)
    log_q = log(q_t_new + 1e-100);
else
    log_q = 0;
end

end