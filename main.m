clear all
close all
clc

 T = 1; %payoff time
 dt = 0.01; %timestep
 n_timesteps = T/ (dt); %number of timesteps
 h = 0.2; %discretization parameter for derivative

 eps = [0.01 0.02 0.04 0.08 0.1 0.15 0.2 0.3 0.4 1]; %model error parameter array
 %model parameters
 lambda = 0.4-log(1./eps);
 sigma = 0.2*ones(1,length(eps));
 sigma_eps = sqrt(eps);

 n_sims = 1000; %number of Monte-Carlo simulations
 MC_iterations = []; %vector of results from Monte Carlo simulations
 G3_estimations = []; %vector of results from estimation formula

 for i = 1 : n_sims
   
    %simulate Compound Poisson process with random jumps
    %first construct a vector of random jumps to apply to the process
    n_jumps = poiss (1,25); %number of jumps taken from a Poisson distribution
    size_jumps = randn(n_jumps,1); %size of jumps taken from a Normal distribution
    jumps = sum(size_jumps); % sum of all jumps to construct 3rd term of simulation process
    
    %simulate the stochastic process
    X_T = lambda.*T + (sigma+sigma_eps).*Wiener(n_timesteps,dt) + jumps.*ones(1,length(eps)); %compute stochastic step
    %calculate payoff
    payoff = payoff_g(X_T);
    %calculate 3rd derivative of payoff
    g3 = (1/(8*h^3))*(payoff_g(X_T(end)-3*h)-8*payoff_g(X_T(end)-2*h)+13*payoff_g(X_T(end)-h)-13*payoff_g(X_T(end)-h)+8*payoff_g(X_T(end)+2*h)-payoff_g(X_T(end)-3*h));

    MC_iterations = [MC_iterations; payoff]; %matrix of Monte Carlo iteration results
    G3_estimations = [G3_estimations; g3]; %vector of results from 3rd derivative of payoff
 end

 E = (1/n_sims).*sum(MC_iterations); %expected value of payoffs
 Eg3 = (1/n_sims).*sum(G3_estimations); %expected value of 3rd derivative of payoffs

 %calculate the errors in Monte Carlo estimation
 errors = zeros(1,length(eps)-1);
 for i = 1 : length(eps)-1
    errors(i) = abs(E(1)-E(i+1));
 end

 %calculate theoretical error 
 eps_cont = linspace (0.01,1,1000);
 error_M = (T/6)*((eps_cont.^2)/6)*Eg3;

 