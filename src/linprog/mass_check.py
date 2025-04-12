from facfac import main
import math

results = open("../../Data/results.txt", "a")

lower = 43000
upper = 43631

save_path_partial = "../../Data/low_factorizations2/"

for i in range(lower, upper+1):
    check = main(i, math.ceil(i/3), save=save_path_partial + "factors_" + str(i) + ".txt")
    if check:
        results.write(f"{i}: SUCCESS\n")
    else:
        results.write(f"{i}: FAILURE\n")