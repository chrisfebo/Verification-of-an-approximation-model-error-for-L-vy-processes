function [ X ] = poiss( numel,lambda )
% generates a vector of random numbers from a Poisson distribution

for i=1:numel
    k=1;
    usave=1;
    usave = usave*rand;
    while usave >= exp(-lambda)
        usave = usave*rand;
        k = k+1;
 end
    X(i) = k;
end

end
