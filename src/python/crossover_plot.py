import numpy as np
import matplotlib.pyplot as plt

data = {}

file = "../../Data/crossover_bounds.txt"

with open(file, 'r') as f:
    for line in f:
        key_str, value_str = line.strip().split(": ")
        key = int(key_str)
        value = eval(value_str)
        data[key] = value

Ns = np.array(list(data.keys()))

transformed_bounds = [(lower - N, upper - N) for N, (lower, upper) in data.items()]
transformed_bounds = np.array(transformed_bounds)

# Plotting the transformed bounds
plt.figure(figsize=(10, 6))
plt.plot(Ns, transformed_bounds[:, 1], label='Upper Bound')

# Add a horizontal line at y = 0
plt.axhline(y=0, color='r', linestyle='--', label='y = 0')

plt.xlabel('N')
plt.ylabel('Bound Difference From N')
plt.title('Bound Difference From N')
plt.legend()
plt.savefig("../../Data/crossover_difference.png")
plt.show()
