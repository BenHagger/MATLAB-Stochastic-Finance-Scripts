function [Stock_Prices,Call,Put, Analytical_Call, Analytical_Put] = Black_Scholes_European_price(start_stock_price,start_time, final_time,risk_free_rate,sigma,M,timestep_size, Exercise_price)
% M : number of simulated paths desired. e.g. do you want to flip M = 10
% coins or M = 10000 coins at once to find the probability of getting a
% head (larger M means more accurate pricing, at cost of computation)

% risk free rate usually given by a 3-month T-Bill
% sigma is the constant IMPLIED volatility, extracted from market option
% prices
% start_price is simply the current price of the stock

% use time input to create a time axis, t
t = start_time:timestep_size:final_time;
% Work out length of trajectory from timestep
N = length(t); % N is the number of timesteps in each trajectory

% duplicate 1xN time vector t, M times to create a M x N matrix for use in 
%vectorised operations. 'time' will look like:
%  t11, t21, ... , tN1
%  t12, t22, ... , tN2
%   ......
%  t1M, t2M, ... , tNM
time = repmat(t,M,1);

%create increments of brownian motion by using sqrt(dt)*N(0,1). 
% By shift theorem, BM_increments ~ N(0,dt)
BM_increments=sqrt(timestep_size)*randn(M,N-1); 

% concate arrays with initial position of zero, gives our set of M random
% walks (plot to see what this is). cumsum(-, 2) sums our arrays along the
% row, as to add up the simulated increments in time - using the steps to
% generate a walk
BM_whole_path=[zeros(M,1) cumsum(BM_increments, 2)]; 

%create MxN array of trajectories. We do this by taking the exact solution
% S_T = S_0 * exp{ (mu - 0.5* sigma^2) * t  + sigma*W_T  }
% Which solves the Black-Scholes SDE for a log-normal random walk:
% dS = mu*S*dt + sigma*S*dW_t
% See    https://en.wikipedia.org/wiki/Geometric_Brownian_motion 
Stock_Prices = start_stock_price*exp((risk_free_rate-0.5*sigma.^2).*time + sigma.*BM_whole_path); 

% we take the final value of S, which is the asset price at expiry, and
% calculate the payoff each  of the M trajectories would have earned on on the contract
payoff_call = max( Stock_Prices(:, end) - Exercise_price , 0);
payoff_put = max( Exercise_price - Stock_Prices(:, end) , 0 );

% We then take the ensemble average of the payoffs, and discount this back
% to time t = 0 from t = T
Call = exp(-risk_free_rate.*(final_time)).*mean(payoff_call);
Put = exp(-risk_free_rate*(final_time)).*mean(payoff_put);

%Calculate Analytical Call and Put price using d1, d2
% this is simply the implementation of the Black-Scholes analytical formula
d1 = (log(start_stock_price./Exercise_price) + (risk_free_rate+0.5*sigma^2)*final_time )/(sigma*sqrt(final_time));
d2 = d1 - sigma*sqrt(final_time);
% normcdf is the cumulative normal distribution function
Analytical_Call = start_stock_price .* normcdf(d1) - Exercise_price.*exp(-risk_free_rate.*final_time).*normcdf(d2);
Analytical_Put = Exercise_price.*exp(-risk_free_rate.*final_time).* normcdf(-d2) - start_stock_price.*normcdf(-d1);
end