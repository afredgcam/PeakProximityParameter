function log_qA = log_q_A1(A,peak,data,t_given,timesNow,mu_A_prior,sig_A2_prior,sig2,sig_y2)

% Computes probability of sampling the 'peak'-th amplitude

% get number of peaks
N = length(t_given);
numInds = length(timesNow);

% make joint prior parameters
[~,len_mu] = size(mu_A_prior);
if len_mu > 1
    mu_A = mu_A_prior;
else
    mu_A = mu_A_prior*ones(1,N);
end
Sig_A = sig_A2_prior*eye(N);

% get signal
F = signalComponents(t_given,timesNow,sig2); % numInds x N

% get marginal likelihood variance
SigmaData = F * Sig_A * F' + sig_y2*eye(numInds); % numInds x numInds

% avoid taking an inverse 
rescaleFactor = linsolve(SigmaData,F * Sig_A); % numInds x N

% compute difference between data and expectation
observedVariation = data - mu_A * F'; % 1 x numInds

% compute posterior mean
posterior_mu_A = mu_A + observedVariation * rescaleFactor; % 1 x N

% compute posterior variance
posterior_Sig_A = Sig_A - (F * Sig_A)' * rescaleFactor; % N x N

% get the non-'peak' quantities
Sig_not_i = posterior_Sig_A([1:peak-1,peak+1:N],[1:peak-1,peak+1:N]); % N-1 x N-1
Sig_i_not_i = posterior_Sig_A(peak,[1:peak-1,peak+1:N]); % 1 x N-1
mu_not_i = posterior_mu_A([1:peak-1,peak+1:N]); % 1 x N-1
A_not_i = A([1:peak-1,peak+1:N]); % 1 x N-1

%%% compute the conditional for 'peak'

% avoid taking an inverse...again!
rescaled_i = Sig_not_i \ Sig_i_not_i'; % N-1 x 1

% get conditional mean
mu_i = posterior_mu_A(peak) + (A_not_i - mu_not_i) * rescaled_i; % 1 x 1

% get conditional variance
sig_i2 = posterior_Sig_A(peak,peak) - Sig_i_not_i * rescaled_i; % 1 x 1

% compute log(sample density)
A_vari = A(peak) - mu_i;
const = -log(2*pi*sig_i2) / 2;
log_qA = const - 0.5 * A_vari^2 / sig_i2;

end