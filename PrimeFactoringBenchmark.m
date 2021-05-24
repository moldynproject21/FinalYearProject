x = 19;
for i = 1:6
    primeNumbers = primes(uint64(2^x));
    compositeNumbers = primeNumbers.*primeNumbers(randperm(numel(primeNumbers)));
    factors = zeros(numel(primeNumbers),2);

    tic;
    parfor idx = 1:numel(compositeNumbers)
        factors(idx,:) = factor(compositeNumbers(idx));
    end
    toc
end