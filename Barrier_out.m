function [price] = Barrier_out(s, E, BU, BD, r, T)
% function is designed to handle Up/Down and out options of a Monte Carlo
% Scheme, only reqiuring an array of simulated asset PATHS, s

% up and out contract
if isempty(BD) == 1 % if no down barrier is given
    disp('Up and Out selected')
    % find values above the barrier and extract the row (pathway)
    [row, ~] = find(s>BU);
    % set any rows with values above barrier to zero
    s(row,:) = 0;
    % create payoff
    payoff = max(s(:,end)-E,0);
    % average payoff to calculate price
    price = exp(-r*T)*mean(payoff);

elseif isempty(BU) == 1 % no up barrier given
    disp('Down and Out selected')
    % same process as before
    [row, ~] = find(s<BD);
    s(row,:) = 0;
    payoff = max(s(:,end)-E,0);
    price = exp(-r*T)*mean(payoff);

else % both barriers given
    disp('Double knock out selected')
    [row,~] = find(s<BD|s>BU);
    s(row,:) = 0;
    payoff = max(s(:,end)-E,0);
    price = exp(-r*T)*mean(payoff);
end
end
