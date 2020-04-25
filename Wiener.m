function [ W ] = Wiener( n,dt )
% generates a random output from a Wiener process of n steps with timestep dt
W=0;
for i = 1:n
dW = sqrt(dt)*randn(1);  % Wiener increments
W = W+dW;          % Brownian path
end
