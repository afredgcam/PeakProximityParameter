function [overlap,underlap] = averageImperfections(meanSamples,ampSamples,interval,sig2)

% Finds the overlap and cancellation (underlap) for each sampled curve; returns mean of each.

% get useful counts
indMax = length(interval);
numSamples = length(meanSamples);

% prepare over- and under-laps for all samples
allOver = zeros([1,numSamples]);
allUnder = zeros([1,numSamples]);

for j = 1:numSamples

    % get quantities for sample j
    t = meanSamples{j};
    A = ampSamples{j};

    % get number of peaks in sample j
    Nj = length(t);
    if Nj == 0
        allOver(j) = 0;
        allUnder(j) = 0;
        continue
    end

    % get total (absolute) peak area
    allAmps = sum(abs(A));
    
    %%%%% overlap %%%%%

    for i1 = 1:Nj-1
        for i2 = i1+1:Nj
            if sign(A(i1)) ~= sign(A(i2))
                % skip if of different signs
                continue
            end
            % label correctly
            if t(i1) < t(i2)
                t_low = t(i1);
                t_high = t(i2);
                A_low = A(i1);
                A_high = A(i2);
            else
                t_low = t(i2);
                t_high = t(i1);
                A_low = A(i2);
                A_high = A(i1);
            end
            % add overlap to total
            xPrime = sig2 * log(abs(A_low/A_high)) / (t_high - t_low) + (t_high + t_low) / 2;
            phi_high = normcdf(xPrime,t_high,sig2^0.5);
            phi_low = normcdf(xPrime,t_low,sig2^0.5,'upper');
            allOver(j) = allOver(j) + abs(A_high)*phi_high + abs(A_low)*phi_low;
        end
    end

    % standardise result
    allOver(j) = allOver(j) / allAmps;

    %%%%% underlap %%%%%

    % get absolute expected curve, |x_k|
    x_k = abs(A * signalComponents(t,interval,sig2)');

    % get approximate area under |x_k|
    total_area = sum(x_k) - (x_k(1) + x_k(indMax)) / 2;

    % get ratio
    allUnder(j) = 1 - (total_area / allAmps);

    %%%%%%%%%%%%%%%%%%%%

end

% get over- and under- laps
overlap = mean(allOver);
underlap = mean(allUnder);

end