function [t_i,A] = neutralStep(t_given,A_given,data,indMax,parameters)

% Proposes no change in peak number; returns the accepted sample

%%% get known parameters

% peak variance
sig2 = parameters.sig2;
% measurement error variance
sig_y2 = parameters.sig_y2;
% prior amplitude mean (same for all peaks)
mu_A_prior = parameters.mu_A_prior;
% prior amplitude variance (same for all peaks)
sig_A2_prior = parameters.sig_A2_prior;
% peak proximity parameter
eps = parameters.eps;


%%% do the fixed-dimension step

% get number of peaks
N = length(t_given);
% vectorise prior mean amplitudes
mu_A_prior_N = mu_A_prior * ones([1,N]);
% start evolving quantities
A = A_given;
t_i = t_given;

for i = 1:N
    % propose sample for t_i and A
    [try_t_i,log_qt] = sampleTimesApart(data,t_i,1:indMax,A,i,eps,sig2,true);
    [try_A,log_qA_new] = sample1Amp_with_prob(data,try_t_i,A,i,1:indMax,mu_A_prior_N,sig_A2_prior,sig2,sig_y2);
    log_qA_old = log_q_A1(A,i,data,t_i,1:indMax,mu_A_prior_N,sig_A2_prior,sig2,sig_y2);
    % get log(acceptance probability)
    log_lik_ratio = logLikDifAll(data,try_t_i,t_i,try_A,A,1:indMax,sig2,sig_y2);
    log_A_ratio = - (try_A(i) - A(i)) * (try_A(i) + A(i) - 2*mu_A_prior) / (2*sig_A2_prior);
    log_accept = log_lik_ratio + log_A_ratio + log_qA_old - log_qA_new + log_qt;

    if log(rand(1)) < log_accept
        t_i = try_t_i; % accept proposed sample
        A = try_A;
    end

end

end


