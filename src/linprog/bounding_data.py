from find_crossover import log_random
from facfac import sieve, lpfac, greedy
import math
import matplotlib.pyplot as plt
import numpy as np

LOW = 10**2 # low end of intitial range
HIGH = 10**5 # high end of initial range
SAMPLE_COUNT = 300 # number of samples to take

def get_bounds(N, iteration_limit=5):
    # get initial upper and lower bounds on T(N)
    Tlb = math.ceil(2*N//7)
    if N == 56:
        Tlb -= 1
    Tub = math.ceil(N/math.exp(1))

    iterations = 0
    while Tub-Tlb > 1 and iterations < iteration_limit:
        T = (Tlb+Tub)//2
        print(T, Tlb, Tub)

        prob = sieve(N, T)

        fact, upper_bound = lpfac(prob)
        if N - upper_bound > 0.001: # check to make sure the upper bound for N is a decent bit below N
            Tub = T
            print("Shifting upper bound")
        else:
            fNT = greedy(prob, fact) # attempt to factorize
            nfactors = sum(fNT.f.values())
            if nfactors >= N:
                Tlb = T
                print("Shifting lower bound")
            else: # cannot conclusively establish a lower bound or an upper bound on this value
                break
        
        iterations += 1
    print(T, Tlb, Tub)

    return Tlb, Tub

samples = [log_random(LOW, HIGH) for i in range(SAMPLE_COUNT)] # sample points from log uniform distribution
samples.append(LOW) # add lower and upper points to the distribution
samples.append(HIGH)
samples.sort()

ITERATIONS = 5
bounds = [get_bounds(sample, iteration_limit=ITERATIONS) for sample in samples]

samples = np.array(samples)
bounds = [(bounds[i][0]/samples[i], bounds[i][1]/samples[i]) for i in range(len(samples))]
bounds = np.array(bounds)

plt.figure(figsize=(15, 5))
plt.scatter(samples, bounds[:, 1], label='Upper Bound')
plt.scatter(samples, bounds[:, 0], label='Lower Bound')
plt.xscale('log')
plt.xlabel('Number of Samples (log scale)')
plt.ylabel('Normalized Bound')
plt.title('Sample Size vs. Normalized Bounds')
plt.legend()
plt.tight_layout()
plt.savefig("output.png")
plt.show()