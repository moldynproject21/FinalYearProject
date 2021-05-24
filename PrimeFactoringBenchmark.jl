using Primes
using Random

### Number Of Primes
x = 19

primeNumbers = Primes.primes((2^x))
compositeNumbers = primeNumbers.*primeNumbers[randperm(length(primeNumbers))]
factors = zeros(length(primeNumbers),2)


t1 = time()
for idx = 1:length(compositeNumbers)
    factors[idx,:] = factor(Vector, compositeNumbers[idx])
end
t2 = time()

println("Seconds Taken: ", (t2-t1))