# plots of t(N)/N with some primes restricted
import matplotlib.pyplot as plt
import math

base = []
values=[]

def read_file(name):
    # Open the file for reading
    with open(name, "r") as file:
        # Read the file line by line
        for line in file:
            # Split the line into parts (assuming numbers are separated by whitespace)
            parts = line.split()
            # Ensure there are at least two numbers on the line
            if len(parts) >= 2:
                base.append(int(parts[len(parts)-2]))
                values.append(int(parts[len(parts)-1]))

def plot2():
    read_file("..\\..\\..\\data\\rearrange\\up_to_2.txt")

    exact = [values[i]/base[i] for i in range(len(base))]
    compare = [3/16 for N in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, exact, label='$t_2(N)/N$' )
    plt.plot(base, compare, linestyle="--", label='$3/16$' )
    plt.title('$t_2(N)/N$')
    plt.xlabel('$N$')
    plt.ylim(0.187,0.189)
    plt.legend()
    plt.grid(True)
    plt.show()

def plot3():
    read_file("..\\..\\..\\data\\rearrange\\up_to_3.txt")

    exact = [values[i]/base[i] for i in range(len(base))]
    compare = [1/4 for N in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, exact, label='$t_{2,3}(N)/N$' )
    plt.plot(base, compare, linestyle="--", label='$1/4$' )
    plt.title('$t_{2,3}(N)/N$')
    plt.xlabel('$N$')
    plt.ylim(0.245,0.255)
    plt.legend()
    plt.grid(True)
    plt.show()


def plot5():
    read_file("..\\..\\..\\data\\rearrange\\up_to_5.txt")

    exact = [values[i]/base[i] for i in range(len(base))]
    compare = [math.ceil(2*N//7)/N for N in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, exact, label='$t_{2,3,5}(N)/N$' )
    plt.plot(base, compare, label='$2/7$', color='brown' )
    plt.title('$t_{2,3,5}(N)/N$')
    plt.xlabel('$N$')
    plt.ylim(0.275,0.286)
    plt.legend()
    plt.grid(True)
    plt.show()

def plot7():
    read_file("..\\..\\..\\data\\rearrange\\up_to_7.txt")

    exact = [values[i]/base[i] for i in range(len(base))]
    compare = [math.ceil(2*N//7)/N for N in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, compare, label='$2/7$', color='brown' )
    plt.plot(base, exact, label='$t_{2,3,5,7}(N)/N$' )
    plt.title('$t_{2,3,5,7}(N)/N$')
    plt.xlabel('$N$')
    plt.ylim(0.285,0.293)
    plt.legend()
    plt.grid(True)
    plt.show()

plot7()