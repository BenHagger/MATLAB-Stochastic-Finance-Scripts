function [Call, S, V] = Heston_Down_Out_exit_prob(k, theta, v0, p, sigma, T, s0, dt_or_N, r, M, E, D)
% k: mean reversion speed of v_t
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
% D: lower barrier
% condition which of the inputs is the number of steps or length of
% timestep
if dt_or_N>1
    N = dt_or_N;
    t = linspace(0,T, N);
    dt = T/N;
else
    dt = dt_or_N;
    t = 0:dt:T; % time axis
    N = length(t); % length of each simulated pathway
end
tic
% function uses Milstein discretization on the CIR model to compute values
% of v_t with the truncation scheme v_t = max{v_t,0} to avoid complex prices
% It deploys the Milstein scheme for s_t at each timestep

V = zeros(M,N);
V(:,1) = v0;
S = zeros(M,N);
S(:,1) = s0;
%precompute random variables for for loop
dZs = sqrt(dt).*randn(M,N-1);
dZv = p*dZs +sqrt(1- p^2).*sqrt(dt).*randn(M,N-1);

%create variable which will be used to zero pathways that hit barrier
ActiveM = true(M,1);
for i = 1:N-1
    %extract current number of non-zero rows
    Rows = find(ActiveM);
    %store increments of brownian motion
    dZv_M = dZv(Rows,i);
    dZs_M = dZs(Rows,i);

    % take previous values and set to constants
    v_p = V(Rows,i);
    s_p = S(Rows,i);

    %Milstein scheme with full truncation
    v_next =v_p - k*(max(v_p,0)-theta)*dt + sigma.*sqrt(v_p).*dZv_M + 0.25*sigma.^2 .*(dZv_M.^2 -dt);
    %v_next =v_p +k*(theta - v_p)*dt +sigma.*sqrt(v_p*dt).*dZv_M + 0.25*sigma.^2 .*dt.*(dZv_M.^2 -1);
    s_next = s_p + (r).*s_p.*dt + sqrt(v_p).*s_p.*dZs_M; %+ 0.25*s_p.*v_p .*(dZs_M.^2 - dt);
    
    %remove any prices which hit barrier from next iteration
    ActiveM(Rows) = s_next>D;
    % Calculate exit probabilities only for active rows
    psi = rand(length(Rows), 1);
    Tprob = exp( (-2).*(1./(sigma.^2 * dt .* (sqrt(v_p).*s_p).^2 )).*(D - s_p).*(D - s_next) );
    
    % Update ActiveM based on the comparison between Tprob and psi
    ActiveM(Rows) = psi > Tprob;

    % Take newly generated values and place into variable arrys % use
    % reflection schemes
    V(Rows,i+1) = max(v_next,0); 
    S(Rows,i+1) = s_next; 
end

payoff_call = max( S(:, end) - E , 0);
% We then take the ensemble average of the payoffs, and discount this back
% to time t = 0 from t = T
Call = exp(-r.*(T)).*mean(payoff_call);
toc
end