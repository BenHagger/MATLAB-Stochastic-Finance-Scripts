function [S, V, Call, Put] = Heston_European_Option(k, theta, v0, p, sigma, T, s0, dt, r, M, E)
% k mean reversion speed of v_t
% theta: long run variance of v_t
% v0: initial variance of v_t
% p: correlation of Brownian Motions
% sigma: volatility of v_t
% T: expiry time of option
% s0: initial price
% dt: timestep 1>>>dt
% mu: long run mean, for risk neutral price, set mu = r
% r: risk free interest rate
% M: number of desired simulations
% E: exercise price
t = 0:dt:T; % time axis
N = length(t); % length of each simulated pathway
tic
% function uses Milstein discretization on the CIR model to compute values
% of v_t with Roger Lord, Remmert Koekkoek, and Dick JC Van Dijk. A comparison of
% biased simulation schemes for stochastic volatility models. (2008)'s full
% truncation scheme avoid complex prices
% It deploys the Milstein scheme for s_t at each timestep

if 2 * k * theta > sigma^2 
% if feller condition satisfied then carry out computation as required
    V = zeros(M,N);
    V(:,1) = v0;
    S = zeros(M,N);
    S(:,1) = s0;
    %precompute random variables for efficiency
    dZs = sqrt(dt).*randn(M,N-1);
    dZv = p*dZs +sqrt(1- p^2).*sqrt(dt).*randn(M,N-1);
    
    for i = 1:N-1
        %store increments of brownian motion
        dZv_M = dZv(:,i);
        dZs_M = dZs(:,i);
        % take previous values and set to constants
        v_p = V(:,i);
        s_p = S(:,i);
        %Milstein scheme with full truncation
        v_next =v_p - k*(max(v_p,0)-theta)*dt + sigma.*sqrt((max(v_p,0)).*dZv_M + 0.25*sigma.^2 .*(dZv_M.^2 -dt);
        s_next = s_p + (r).*s_p.*dt + sqrt(v_p).*s_p.*dZs_M + 0.25*s_p.*v_p .*(dZs_M.^2 - dt);
        % Take newly generated values and place into variable arrys % use
        % reflection schemes
        V(:,i+1) = max(v_next,0); 
        S(:,i+1) = s_next; 
    end
    payoff_call = max( S(:, end) - E , 0);
    payoff_put = max( E - S(:, end) , 0 );
    % We then take the ensemble average of the payoffs, and discount this back
    % to time t = 0 from t = T
    Call = exp(-r.*(T)).*mean(payoff_call);
    Put = exp(-r*(T)).*mean(payoff_put);
else 
    disp('Error: Feller Condition not satisfied')
end
toc
end
