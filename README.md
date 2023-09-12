# MATLAB-Stochastic-Finance-Scripts
Here we have:\
**Barrier_out** - a function applied to any generated Monte Carlo prices, which filter out stock prices which have hit the barrier, and subsequently discards the payoff to price down-and-out Eurpoean Call options. 

**Black_Scholes_European_price** - gives the price of a European Call and Put option using vectorised Monte Carlo methods, and checks the result using the Analytical Black-Scholes prices. 

**Heston_European_Option_price** - Uses the Heston stochastic volatility model in a Monte Carlo simulation in conjunction with a fully truncated Milstein discretization to accurately price European options. 

**Heston_Down_out** - Within the Monte Carlo for loop, a condition is implemented to make asset payoffs zero if a barrier is hit for a down-and-out European Call option. As a result, the number of computations in each iteration of the for-loop decreases every time a pathway hits the barrier. 

**Heston_Down_out_exit_prob** - Alongside the barrier condition, the theoretical probability of a continuous time stock exiting a discretely-monitored barrier is calculated, and compared with uniformly randomly generated variables. If the stock is deemed to have breached the barrier between timesteps, the payoff is considered zero and the number of computations decreases.
