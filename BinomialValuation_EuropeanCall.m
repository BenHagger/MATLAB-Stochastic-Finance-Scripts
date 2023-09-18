function [call_loop, call_vector, call_series] = BinomialValuation_EuropeanCall(stockPrice, strikePrice, maturity, interestRate, volatility, timeStep)

    % Binomial valuation of a European call option. Once optimised, a European
    % put can be calculated, and the two can be altered into american options
    % by implementing the condition

    % The magnitudes 'up' and 'down' are from 'An introduction to Financial
    % Option Valuation', Higham (2004), p.154, which correspond to expecting
    % the stock to move one standard deviation, sigma*sqrt(dt), away from the mean at that
    % timestep, in either the lower or upper direction


    % Start timer for loop method - unvectorised loops
    tic;

    % Calculated parameters
    num_steps = maturity/timeStep;
    up_mag = exp(volatility*sqrt(timeStep) + (interestRate-0.5*volatility^2)*timeStep);
    down_mag = exp(-volatility*sqrt(timeStep) + (interestRate-0.5*volatility^2)*timeStep); 
    discountFactor = exp(-interestRate * timeStep); 
    riskNeutralProb = (1/discountFactor - down_mag)/(up_mag-down_mag);  

    % Initialize array for option values
    PriceOption = zeros(num_steps+1, num_steps+1);

    % Calculate option values at maturity using loop
    for k = 1:(num_steps+1)
        % initialise the final payoffs of the options
        PriceOption(num_steps+1,k) = max(stockPrice*up_mag^(k-1)*down_mag^(num_steps-k+1) - strikePrice, 0); 
    end
    for m = num_steps:-1:1  % loop through each point in time - represented by the columns
        for k = 1:m % loop through each branch at a specific point in time - each of the potential movements of the stock from the previous value
            % loop through 
            PriceOption(m,k) = discountFactor*(riskNeutralProb*PriceOption(m+1,k+1) + (1-riskNeutralProb)*PriceOption(m+1,k));
        end
    end
    call_loop = PriceOption(1,1);

    % End timer for loop method
    time_loop = toc

    % Start timer for vectorized method
    tic;

    % Initialize array for option values
    optionValues = zeros(num_steps+1, num_steps+1);

    % Calculate option values at maturity using vectorized operations
    stepVector = 0:num_steps;
    optionValues(end, :) = max(stockPrice * up_mag.^stepVector .* down_mag.^(num_steps - stepVector) - strikePrice, 0);

    % Backward induction to calculate option values, but with inner loop
    % vectorised
    for m = num_steps:-1:1
        optionValues(m,1:m) = discountFactor*(riskNeutralProb*optionValues(m+1, 2:m+1)+ (1-riskNeutralProb)*optionValues(m + 1, 1:m));
    end

    % Final option value for vectorized method
    call_vector = optionValues(1, 1);

    % End timer for vectorized method
    time_vector = toc

    %%%  series approximation method - using the power series approximation
    %%%  of the binomial formula
    %%% This method is ONLY VALID FOR MATURITY<=1. Any higher maturity
    %%% produces infinite values due to the size of the factorial
    %%% increasing exponentially when the number of steps increases.
    %%% min size dt is 0.001, any smaller and call price is NaN
    tic;
    % Stirling's approximation for factorial using an anonymous function
    stirling = @(x) x*log(x) + (1-x)*log(1/(1-x));

    % initialise vector to store sum
    B = zeros(1, num_steps);
    nCk_approx = zeros(1, num_steps);
    
    % loop over all paths
    for i = 1:num_steps
        % Stirling's approximation for nchoosek
        warning off
        nCk_approx(i) = exp(num_steps*stirling(i/num_steps));
        B(i) = nchoosek(num_steps,i) .* (riskNeutralProb)^i .* (1-riskNeutralProb)^(num_steps - i) .* max(up_mag^i .* down_mag^(num_steps - i) .* stockPrice - strikePrice, 0);
    end
    warning on
    call_series = exp(-interestRate*maturity)*sum(B);

    % End timer for series approximation method
    time_series = toc

end
