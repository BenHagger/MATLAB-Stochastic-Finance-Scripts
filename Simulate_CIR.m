function r_t = simulate_CIR(alpha, sigma, beta, r0, dt, t0, tend, M)
% beta: Mean interest rate
% r0: Initial position
% dt: Timestep
% t: Time vector
% N: Number of iterations per pathway
% M: Number of simulated pathways

%time axis based on timestep
t = t0:dt:tend;
numSigma = length(sigma);

% Initialize arrays
is_fuller_satisfied = 2 * alpha * beta > sigma.^2;
%3d array to store multiple sigma - this way we can perform vectorised
%operations on multiple values of sigma without the need for a for loop.
%Making the process faster and more efficient
r_t = zeros(M, N, numSigma);
r_t(:, 1, :) = r0;

%reshape sigma to 3d array for vector operstions
sigma_r = reshape(sigma, [1, 1, numSigma]);
for i = 1:N-1 % Loop through r_t
    % creates increment of BM for each trajectory and for each sigma
    dZ_t = sqrt(dt) .* randn(M, numSigma);
    %take previous value of r_t and store in a temporary vector
    r_p = r_t(:, i, :);
    % again reshape dZ_t into 3d matrix
    dZ_t_r = reshape(dZ_t, [M, 1, numSigma]);
    % Milstein Discretization for next value of r_t
    r_next = r_p+ alpha.*(beta -r_p).*dt+sigma_r.*sqrt(r_p).*dZ_t_r + 0.25.*sigma_r.^2 .*(dZ_t_r.^2 -dt);
    %implement condition to ensure no comples values of r_t
    r_t(:, i+1, :) = max(r_next, 0);
end
end