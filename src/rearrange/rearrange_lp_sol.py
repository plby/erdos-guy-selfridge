import math
from sympy import primerange
import sys

def main(P: int, N: int, T: int):
    primes = list(primerange(2, P+1))

    divided_factors = []
    for i in range(1, N+1):
        t = i
        # Divide out the small primes
        for p in primes:
            while (t % p) == 0:
                t //= p
        divided_factors.append(t)

    divided_factors.sort()

    # Verify header of lp_solve solution
    line = sys.stdin.readline()
    assert line == "\n", f"Unexpected line (expected empty line): {line!r}"
    line = sys.stdin.readline()
    assert line.startswith("Value of objective function: "), f"Unexpected line (expected \"Value ...\"): {line!r}"
    line = sys.stdin.readline()
    assert line == "\n", f"Unexpected line (expected empty line): {line!r}"    
    line = sys.stdin.readline()
    assert line == "Actual values of the variables:\n", f"Unexpected line (expected \"Actual values ...\"): {line!r}"

    # Reconstruct factors to multiply by
    multiply_factors = []
    for line in sys.stdin:
        variable, value = line.split()
        value = int(value)
        assert variable.startswith("a")
        n = int(variable[1:])
        t = n
        for p in primes:
            while (t % p) == 0:
                t //= p
        assert t == 1, f"Variable {variable} which corresponds to n={n} was not {P}-smooth!"
        for i in range(value):
            multiply_factors.append(n)

    assert len(divided_factors) == len(multiply_factors), f"Unexpected number of factors."
    factors = [divided_factor * multiply_factor for divided_factor, multiply_factor in zip(divided_factors, multiply_factors)]
    print(N, N, T)
    for factor in factors:
        print(factor)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
            prog='rearrange_lp_sol',
            description='Rearrange small prime factors from the lp_solve solution of a rearrange_lp problem.')
    parser.add_argument('P', type=int)
    parser.add_argument('N', type=int)
    parser.add_argument('T', type=int, nargs='?', default=0)

    args = parser.parse_args()

    P = args.P
    N = args.N
    T = args.T or math.ceil(N/4)

    main(P, N, T)
