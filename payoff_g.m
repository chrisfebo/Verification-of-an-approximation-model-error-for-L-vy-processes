function [ X ] = payoff_g( x )
% this function calculates the payoff of a digital option with strike at 0

X = zeros(1,length(x));

for i = 1 : length(x)
  if x(i) > 0
    X(i) = 1;
end

end

end
