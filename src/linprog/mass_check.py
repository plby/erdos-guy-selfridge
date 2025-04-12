from facfac import main
import math

lower = 38000
upper = 40000

results = open(f"../../Data/conjecture2_successes/results_{lower}-{upper}.txt", "a")
save_path_partial = "../../Data/conjecture2_factorizations/"

for i in range(lower, upper):
    check = main(i, math.ceil(i/3), save=save_path_partial + f"factors_{i}.txt")
    if check:
        results.write(f"{i}: SUCCESS\n")
    else:
        results.write(f"{i}: FAILURE\n")