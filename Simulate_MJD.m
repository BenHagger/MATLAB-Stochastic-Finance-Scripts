function [t, S] = Simulate_MJD(lambda, r, sigmaBS, muJD, sigmaJD, S_0, dt, t_0, t_end, nsim)
% sigmaJD is the variance of the jump magnitude
% muJD is the mean jump size
% lambda is the intensity i.e. the mean number of jumps per unit time

% this script uses the definition of the price of a lognormally distributed
% stock price experiencing jumps as found in Merton (1976)
% S_t = S_0 exp{ (mu_BS - lambda*mu_JD - 0.5*sigmaBS^2)*t + sigmaBS*W_t + sum_{i=1}^{i = N_t} U_i  }
% Where N_t is a compount poisson process denoting the random number of
% jumps in our time interval
% U_i is the magnitude of each jump

t = t_0:dt:t_end;
N = length(t);

% Merton (1976)
% create npaths x N array of jump locations using compound poisson process, and occasionally
% will have two jumps instead of one. We can either keep this the same and use it to represent a big
% jump, or we can turn values greater than one into 1 as Merton states
% p(dN_t>1)~= 0 
Jump_location = poissrnd(lambda*dt, [nsim, N-1]);
% find any jump location where the number of jumps is larger than 1
rare_jumps = find(Jump_location>1);
% convert jump value to one jump only
Jump_location(rare_jumps) = 1;

% as the model calls for the addition of normally distributed *independent*
% jumps, we can simply create an array of N(muJ, sigmaJ^2) jumps and filter
% them using our jump locations to filter the array
% will scale and shift normal distribution so J~N(muJ, sigmaJ^2)
Jump_magnitude = muJD.*ones(nsim,N-1) + sigmaJD*randn(nsim,N-1);

% can now filter jump magnitude values to create our jump process
Jumps = Jump_magnitude.*Jump_location;
% create a pathway of jumps by adding in time, whilst concating with an
% initial zero to ensure S_0 = S_0 *1 remains true
Jump_pathway = [zeros(nsim,1), cumsum(Jumps, 2)];

% to model stock prices must create a pathway of Brownian motion
BM_pathway = [zeros(nsim,1), cumsum(sqrt(dt).*randn(nsim,N-1), 2)];
% using the definition of stock price experiencing jumps
S = S_0.*exp((r - lambda.*muJD - 0.5*sigmaBS^2).*t_end + sigmaBS.*BM_pathway + Jump_pathway);
end







