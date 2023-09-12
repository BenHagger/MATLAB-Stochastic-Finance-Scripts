function [S,Call, Analytical_Call] = Black_Scholes_Down_Out_exit_prob(s0,T,r,sigma,M,dt_or_N, E, D)
% s0 initial price
% Tspan: [t0,T]
% r: risk free deterministic growth
% sigma: volatility of random fluctuations
% M: number of simulations
% dt: timestep, recommend small dt >=0.01 
% S: MxN output of trajectories
% script uses solution of lognormal random walk to directly simulate
% trajectories of Geometric Brownian motion, then calculates the payoff at
% expiry and takes ensemble average to estimate a call for a European put
% and call.
% E: can be a vector as it is not used in path generation
% D: lower barrier which if price goes below payoff is zero
tic
if dt_or_N>1
    N = dt_or_N;
    t = linspace(0,T, N);
    dt = T/N;
else
    dt = dt_or_N;
    t = 0:dt:T; % time axis
    N = length(t); % length of each simulated pathway
end

time = repmat(t,M,1);
%create an array of wiener process, W, using the GPU
dW=sqrt(dt)*randn(M,N-1); % create array of random variables
    % from normal distribution 
% concate arrays with initial position of zero
W=[zeros(M,1) cumsum(dW, 2)]; 
%create MxN array of trajectories
S = s0*exp((r-0.5*sigma.^2).*time + sigma.*W);
%eliminate rows which hit barrier
[row, ~] = find(S<D);
S(row,:) = 0;

% calculate exit probabilities 
psi = rand([ M,N-1]);
Tprob = exp( (-2).*(1./(sigma.^2 * dt .* S(:,1:end-1).^2 )).*(D - S(:,1:end-1)).*(D - S(:,2:end)) );

% eliminate trajectories which have exited
[row2, ~] = find(psi < Tprob);
S(row2, :) = 0;
% we take the final value of S, which is the asset price at expiry, and
% calculate the payoff each trajectory would have earned on on the contract
payoff_call = max( S(:, end) - E , 0);

% We then take the ensemble average of the payoffs, and discount this back
% to time t = 0 from t = T
Call = (exp(-r.*(T)).*mean(payoff_call));

%Calculate Analytical Call and Put price using d1, d2
d1 = (log(s0./E) + (r+0.5*sigma^2)*T )/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
ECall = s0 .* normcdf(d1) - E.*exp(-r.*T).*normcdf(d2);
%Analytical barrier price
d1B = (log((D^2./s0)./E) + (r+0.5*sigma^2)*T )/(sigma*sqrt(T));
d2B = d1B - sigma*sqrt(T);
scaled_call = (s0./D)^(1 - 2*r/(sigma^2)) .* ((D^2./s0) .* normcdf(d1B) - E.*exp(-r.*T).*normcdf(d2B));
Analytical_Call = ECall - scaled_call;
toc
end