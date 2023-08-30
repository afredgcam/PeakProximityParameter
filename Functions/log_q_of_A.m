function log_qA = log_q_of_A(A,signalMeasurements,timesMean,timesNow,mu_A_prior,sig_A2_prior,sig2,sig_y2)

% Computes probability of the samples of all the amplitudes

% get number of peaks
N = length(timesMean);
numInds = length(timesNow);

if N == 0

    % no peaks from which to get log_qA
    log_qA = 0;

else

    % make joint prior parameters
    [~,len_mu] = size(mu_A_prior);
    if len_mu > 1
        mu_A = mu_A_prior;
    else
        mu_A = mu_A_prior*ones(1,N);
    end
    Sig_A = sig_A2_prior*eye(N);
    
    % get signal
    F = signalComponents(timesMean,timesNow,sig2); % numInds x N
    
    % get marginal likelihood variance
    SigmaData = F * Sig_A * F' + sig_y2*eye(numInds); % numInds x numInds
    
    % avoid taking an inverse 
    rescaleFactor = linsolve(SigmaData,F * Sig_A); % numInds x N
    
    % compute difference between data and expectation
    observedVariation = signalMeasurements - mu_A * F'; % 1 x numInds
    
    % compute posterior mean
    posterior_mu_A = mu_A + observedVariation * rescaleFactor; % 1 x N
    
    % compute posterior variance
    posterior_Sig_A = Sig_A - (F * Sig_A)' * rescaleFactor; % N x N
    
    % compute log(sample density)
    A_vari = A - posterior_mu_A;
    const = -N * log(2*pi) / 2 - log(det(posterior_Sig_A)) / 2;
    log_qA = const - 0.5 * A_vari * (posterior_Sig_A \ A_vari'); % A \ b = inv(A) * b

end

end