function [t_i,A] = birthStep(t_given,A_given,data,indMax,maxN,parameters)

% Proposes creation of new peak; returns accepted sample

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
% birth probability
p_birth = parameters.p_birth;
% death probability
p_death = parameters.p_death;
% mean peak number
meanN = parameters.meanN;

% correct birth/death probabilities, if needed
if isempty(t_given)
    p_birth = 1;
elseif length(t_given) == maxN - 1
    p_death = p_death + p_birth;
end



%%%%% do birth-proposal step %%%%%


% augment old sample
t_aug = [t_given,0];
A_aug = [A_given,0];
N = length(t_aug);
% sample new peak's delay, and get log_q including prior ratio term
[try_t_i,log_qt_inc_prior] = sampleTimes_with_Prior(data,t_aug,1:indMax,A_aug,N,eps,sig2,true);
% sample all amplitudes
mu_A_prior_N = mu_A_prior * ones([1,N]);
[try_A,log_qA_birth] = sampleAmps_with_prob(data,try_t_i,1:indMax,mu_A_prior_N,sig_A2_prior,sig2,sig_y2);

%%% compute log(acceptance probability)

% compute log(likelihood ratio)
log_lik_ratio = logLikDifAll(data,try_t_i,t_given,try_A,A_given,1:indMax,sig2,sig_y2);

% compute log(A prior ratio)
try_A_var_ss = dot(try_A - mu_A_prior_N,try_A - mu_A_prior_N); % sum{(Ai-mi)^2} for proposal
A_var_ss = dot(A_given - mu_A_prior_N(1:N-1),A_given - mu_A_prior_N(1:N-1)); % sum{(Ai-mi)^2} for old sample
log_A_ratio = -log(2*pi*sig_A2_prior) / 2 - (try_A_var_ss - A_var_ss) / (2*sig_A2_prior);

% compute log(N prior ratio)
log_N_ratio = log(meanN) - log(N);

% compute log q(u,A|ü,Ä) and log q(ü,Ä|u,A)
log_q_birth = log_qA_birth + log_qt_inc_prior + log(p_birth);
log_qA_death = log_q_of_A(A_given,data,t_given,1:indMax,mu_A_prior_N(1:N-1),sig_A2_prior,sig2,sig_y2);
log_death_prob = peakDeathProb(try_A);
log_death_prob = log_death_prob(N);
log_q_death = log_qA_death + log_death_prob + log(p_death);

% combine all of the probabilities
log_prob_accept = log_lik_ratio + log_A_ratio + log_N_ratio + log_q_death - log_q_birth;

%%% accept or reject sample

if log(rand(1)) < log_prob_accept
    % accept
    t_i = try_t_i;
    A = try_A;
else
    % reject
    t_i = t_given;
    A = A_given;
end

end