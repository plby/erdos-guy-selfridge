# code to find bounds on the largest factorization with all factors >= N/3 in a particular region

from smoothfac import main
from math import ceil
from typing import Optional

LOW = 40000 # low end of range
HIGH = 45000 # high end of range
STEP = 1 # steps to skip

results = open(f"output.txt", "w")

# loop through all N in the range
for N in range(LOW, HIGH+1, STEP):
    lower_bound, upper_bound = main(N, ceil(N/3), None, None) # bound number of factors
    results.write(f"{N}: [{lower_bound}, {upper_bound}]\n")
    if N % 100 == 0: # occasionally flush the output so that the file is updated
        results.flush()