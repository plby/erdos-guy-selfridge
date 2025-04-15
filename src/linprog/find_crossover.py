from facfac import main
import numpy as np
import math

output = open(f"crossover_output.txt", "w")

# generate a random number from the log uniform distribution between two bounds
def log_random(lower, upper):
    lower_log = math.log(lower)
    upper_log = math.log(upper)
    random_log = np.random.uniform(lower_log, upper_log)
    random = math.exp(random_log)
    return math.floor(random) # convert to integer

LOW = 10**1 # low end of intitial range
HIGH = 10**6 # high end of initial range
SAMPLE_COUNT = 50 # number of samples to take
CUTOFF = 25 # number of successes to cut off after

def f(x):
    return main(x, math.ceil(x/3)) # change T here

def last_failure(lower, upper):
    output.write(f"Searching range [{lower}, {upper}]\n")

    # if the range is small enough to just examine the entire range, do that
    if upper - lower + 1 <= SAMPLE_COUNT:
        highest_failure = lower
        for sample in range(lower, upper+1):
            success = f(sample)
            if not success:
                highest_failure = sample
        return highest_failure
    else:
        samples = [log_random(lower, upper) for i in range(SAMPLE_COUNT)] # sample points from log uniform distribution
        samples.append(lower) # add lower and upper points to the distribution
        samples.append(upper)
        samples.sort()

        # cycle through samples, testing each point. for each failure, adjust the lower bound of the next range to be at that point
        # for each success, wait until there have been CUTOFF successes in a row, and then halt, taking the last point as the next upper bound
        highest_failure = lower
        contiguous_successes = 0
        results = []
        for sample in samples:
            success = f(sample)
            if success:
                contiguous_successes += 1
            else:
                highest_failure = sample
                contiguous_successes = 0
            results.append(success)
            if contiguous_successes >= CUTOFF:
                break
        output.write(str([(samples[i], results[i]) for i in range(len(results))]) + "\n")
        output.flush()

        # recurse on the restricted range
        return last_failure(highest_failure, sample)

print(last_failure(LOW, HIGH))