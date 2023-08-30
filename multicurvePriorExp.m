%%%%% This is for comparing  RJ-MCMC %%%%%

% NOTE: amplitudes are coeff's of Gaussian peaks, not their max height!!

clear
% may need '\' instead of '/' if on windows/linux
addpath './Functions'; addpath './Plots'; addpath './Synthetic Pulse Data';
rng(1)

%%%%% changable parameters %%%%%

% get interval
indMax = 200;
interval = 1:indMax;

% set the design choices
designChoices.numSamples = 300; % number of samples
designChoices.burnin = 100; % burn-in period
samplesUsed = designChoices.numSamples - designChoices.burnin;
designChoices.indMax = indMax; % data length
designChoices.maxN = 1000; % maximum number of peaks

% set the parameters
sig2_tru = 20;
parameters.sig2 = sig2_tru; % peak width
parameters.sig_y2 = 1; % measurement noise variance
parameters.mu_A_prior = 0; % prior amplitude mean (for each peak)
parameters.sig_A2_prior = 100 * 2*pi*sig2_tru; % prior amplitude variance
parameters.p_birth = 0.1; % birth step probability
parameters.p_death = 0.1; % death step probability
parameters.meanN = 5; % (truncated) poisson mean for peak number

% save true parameters separately
tru_parameters = parameters;

% mismatch peak variance, if desired
sig2 = sig2_tru; % set to sig2_tru if want to match
parameters.sig2 = sig2;
parameters.eps = 0.5 * sig2^0.5; % peak proximity parameter

% choose number of data sets
numSets = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare metrics
metrics = zeros([5,2+1,numSets]); % numMetrics x (numMethods + 1) x numCurves <–– +1 for ground truth

% prepare progress bar
progress = waitbar(0, 'Starting');

for k = 1:numSets

    % make data to use
    [y_k,x_k,t_tru,A_tru] = makeSyntheticData(['used_set_',num2str(k)],interval,tru_parameters,designChoices.maxN);
    
    % run RJ-MCMC method with proposed prior
    tic
    [meansMCMC,ampsMCMC] = runRJMCMC(y_k,designChoices,parameters,false); % uses eps to stop close peaks
    timeTaken = toc;
    avgCurveMCMC = averageCurve(meansMCMC,ampsMCMC,interval,sig2);
    rssMCMC = (sum((avgCurveMCMC - x_k).^2,"all"))^0.5;
    [overMCMC,underMCMC] = averageImperfections(meansMCMC,ampsMCMC,interval,sig2);
    N_MCMC = zeros([1,samplesUsed]);
    for j = 1:samplesUsed
        N_MCMC(j) = length(meansMCMC{j});
    end
    metrics(:,1,k) = [rssMCMC,overMCMC,underMCMC,mean(N_MCMC),timeTaken];
    
    % run RJ-MCMC method without proprosed prior
    tic
    [meansClose,ampsClose] = closeRJMCMC(y_k,designChoices,parameters,false); % eps unused
    timeTaken = toc;
    avgCurveClose = averageCurve(meansClose,ampsClose,interval,sig2);
    rssClose = (sum((avgCurveClose - x_k).^2,"all"))^0.5;
    [overClose,underClose] = averageImperfections(meansClose,ampsClose,interval,sig2);
    N_Close = zeros([1,samplesUsed]);
    for j = 1:samplesUsed
        N_Close(j) = length(meansClose{j});
    end
    metrics(:,2,k) = [rssClose,overClose,underClose,mean(N_Close),timeTaken];
    
    % get metrics for ground truth
    [overTru,underTru] = averageImperfections({t_tru},{A_tru},interval,sig2);
    metrics(:,3,k) = [0,overTru,underTru,length(t_tru),0];

    % update progress bar
    donePercent = floor(100 * k / numSets);
    waitbar(donePercent/100,progress,sprintf('Progress: %d %%',donePercent));

end

% remove complete progress bar
delete(progress)

%%% save the metrics %%%
avgMetrics = mean(metrics,3)';
save('Plots/priorExpMetrics.mat',"avgMetrics")
disp(avgMetrics)


% use some actually decent colours!!
colours = {'#484D7C','#D897F8','#C6FFEA','#FFC2FD','#EEEE65','#00B9B2','#9E90DB','#4600AE'};
numColours = length(colours);
green = '#00F090';
true_green = '#005030';
red = '#F00090';
pink = '#E8A8D2';
purple = '#4600AE';
blue = '#4080FF';
J_use = length(meansMCMC);

%%% make desired plots for final set

% neatly plot estimate with echo curve and ground truth

% set figure and axes
fig1 = figure;
ax1 = axes();
ax1.FontSize = 16;

hold on
% make the plotting handles
h1(1) = plot(interval,y_k,'DisplayName','data',LineWidth=3,Color=green,LineStyle=':');
% h1(2) = plot(interval,x_k,'DisplayName','true curve',LineWidth=3,Color=blue,LineStyle=':');
h1(2) = plot(interval,avgCurveMCMC,'-','Color',red,'DisplayName','\delta used','LineWidth',2);
h1(3) = plot(interval,avgCurveClose,'-','Color',purple,'DisplayName','no \delta','LineWidth',2);
legend(h1,'FontSize',24,'Location','best');

% set axis scales and labels
xlabel('time index')
ylabel('return')
hold off

% save the plot and remove it from screen
saveas(ax1,'Plots/priorExp1.fig')
% delete(fig1)


% neatly plot histogram of sampled N with ground truth

% set figures and axes
fig2 = figure;
ax2 = axes();
ax2.FontSize = 16;

hold on
% make the plotting handles
h2(1) = xline(length(A_tru),'LineStyle','--','Color',green,'LineWidth',3,'DisplayName','true count');
h2(2) = histogram(N_MCMC,'FaceColor',red,'DisplayName','\delta used');
h2(3) = histogram(N_Close,'FaceColor',purple,'DisplayName','no \delta');
legend(h2,'FontSize',24,'Location','best');

% set axis scales and labels
xlabel('number of peaks');
ylabel('frequency');
hold off

% save the plot and remove it from screen
saveas(ax2,'Plots/priorExp2.fig')
% delete(fig2)


