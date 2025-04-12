from facfac import main
import math

lower = 42000
upper = 42999

results = open(f"../../Data/results_{lower}-{upper}.txt", "a")
save_path_partial = "../../Data/low_factorizations/"

for i in range(lower, upper+1):
    check = main(i, math.ceil(i/3), save=save_path_partial + f"factors_{i}.txt")
    if check:
        results.write(f"{i}: SUCCESS\n")
    else:
        results.write(f"{i}: FAILURE\n")