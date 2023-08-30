function [y,x,t,A] = makeSyntheticData(filename,interval,parameters,maxN)

% Creates a synthetic data set and saves it in the relevant folder

% extract necessary parameters
sig2 = parameters.sig2; % peak variance
sig_y2 = parameters.sig_y2; % measurement noise variance
meanN = parameters.meanN; % (poisson) mean number of peaks
mu_A = parameters.mu_A_prior; % peak mean amplitude
sig_A2 = parameters.sig_A2_prior; % peak amplitude variance
K = length(interval); % number of data points

% sample number of peaks
N = min([poissrnd(meanN),maxN]);

% sample peak locations
t = rand([1,N]) * (interval(K) - interval(1)) + interval(1);

% sample peak amplitudes
A = mu_A + randn([1,N]) * sig_A2^0.5;

% make sure no peak is 'too small'
% A = A + sig_A2^0.5 * sign(A);

% get true curve
x = A * signalComponents(t,interval,sig2)';

% sample data
y = x + randn([1,K]) * sig_y2^0.5;

% save the data set
allData = {y,t,A};
save(['Synthetic Pulse Data/',filename],"allData")

end