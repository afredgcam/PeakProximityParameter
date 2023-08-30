function [timesMean,amplitudes] = closeRJMCMC(data,designChoices,parameters,wantingDiagnostics)

% Gets all samples, discarding the burnt in ones

% separate design choices
indMax = designChoices.indMax;
numSamples = designChoices.numSamples;
burnin = designChoices.burnin;
maxN = designChoices.maxN;

% get step-type probabilities
p_birth = parameters.p_birth;
p_death = parameters.p_death;

% prepare output objects
amplitudes = cell([1,numSamples-burnin]);
timesMean = cell([1,numSamples-burnin]);
        
% initialise with zero peaks
t_i = zeros([1,0]);
A = zeros([1,0]);

if wantingDiagnostics
    %%% set some diagnostic objects, to be displayed
    births_proposed = 0;
    deaths_proposed = 0;
    births_accepted = 0;
    deaths_accepted = 0;
end


%%%%% run the RJ-MCMC method %%%%%


for j = 1:numSamples % loop over samples
    
    % sample for step-type
    u = rand();
    % select corresponding step
    if ((u < p_birth) || isempty(A)) && (length(A) < maxN)
        if wantingDiagnostics
            births_proposed = births_proposed + 1;
            N = length(A);
        end
        % do a birth step
        [t_i,A] = birthClose(t_i,A,data,indMax,maxN,parameters);
        if wantingDiagnostics
            if length(A) > N
                births_accepted = births_accepted + 1;
            end
        end
    elseif u - p_birth < p_death
        if wantingDiagnostics
            deaths_proposed = deaths_proposed + 1;
            N = length(A);
        end
        % do a death step
        [t_i,A] = deathClose(t_i,A,data,indMax,maxN,parameters);
        if wantingDiagnostics
            if length(A) < N
                deaths_accepted = deaths_accepted + 1;
            end
        end
    else
        % do a fixed-dimension step
        [t_i,A] = neutralClose(t_i,A,data,indMax,parameters); % neutralStep_old one may be getting better results...
    end

    if j > burnin
        % if past burn-in period, store this sample
        amplitudes{j-burnin} = A;
        timesMean{j-burnin} = t_i;
    end
    
end

if wantingDiagnostics
    disp(['births proposed ',num2str(births_proposed),' - ',num2str(births_accepted),' births accepted'])
    disp(['deaths proposed ',num2str(deaths_proposed),' - ',num2str(deaths_accepted),' deaths accepted'])
end

end