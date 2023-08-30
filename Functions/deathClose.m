function [t_i,A] = deathClose(t_given,A_given,data,indMax,maxN,parameters)

% Proposes death of a peak; returns accepted sample

%%% get known parameters

% peak variance
sig2 = parameters.sig2;
% measurement error variance
sig_y2 = parameters.sig_y2;
% prior amplitude mean (same for all peaks)
mu_A_prior = parameters.mu_A_prior;
% prior amplitude variance (same for all peaks)
sig_A2_prior = parameters.sig_A2_prior;
% birth probability
p_birth = parameters.p_birth;
% death probability
p_death = parameters.p_death;
% mean peak number
meanN = parameters.meanN;

% correct birth/death probabilities, if needed
if length(t_given) == 1
    p_birth = 1;
elseif length(t_given) == maxN
    p_death = p_death + p_birth;
end



%%%%% do death-proposal step %%%%%


% get current number of peaks 
N = length(t_given);
mu_A_prior_N1 = mu_A_prior * ones([1,N]);
mu_A_prior_N2 = mu_A_prior * ones([1,N-1]);
% sample the peak to delete
[peak,log_death_prob] = sampleDeathPeak(A_given);
% get remaining peaks
try_t_i = t_given(1:N ~= peak);
[try_A,log_qA_death] = sampleAmps_with_prob(data,try_t_i,1:indMax,mu_A_prior_N2,sig_A2_prior,sig2,sig_y2);

%%% compute log(acceptance probability)

% compute log(likelihood ratio)
log_lik_ratio = logLikDifAll(data,try_t_i,t_given,try_A,A_given,1:indMax,sig2,sig_y2);

% compute log(A prior ratio)
try_A_var_ss = dot(try_A - mu_A_prior_N2,try_A - mu_A_prior_N2); % sum{(Ai-mi)^2} for proposal
A_var_ss = dot(A_given - mu_A_prior_N1,A_given - mu_A_prior_N1); % sum{(Ai-mi)^2} for old sample
log_A_ratio = log(2*pi*sig_A2_prior) / 2 - (try_A_var_ss - A_var_ss) / (2*sig_A2_prior);

% compute log(N prior ratio)
log_N_ratio = log(N) - log(meanN);

% compute log q(u,A|ü,Ä) and log q(ü,Ä|u,A)
log_q_death = log_qA_death + log_death_prob + log(p_death);
log_qA_birth = log_q_of_A(A_given,data,t_given,1:indMax,mu_A_prior_N1,sig_A2_prior,sig2,sig_y2);
try_A_aug = [try_A(1:peak-1),0,try_A(peak:N-1)];
log_qt_inc_prior = log_q_anyT(data,t_given,1:indMax,try_A_aug,peak,sig2);
log_q_birth = log_qA_birth + log_qt_inc_prior + log(p_birth);

% combine all of the probabilities
log_prob_accept = log_lik_ratio + log_A_ratio + log_N_ratio + log_q_birth - log_q_death;

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