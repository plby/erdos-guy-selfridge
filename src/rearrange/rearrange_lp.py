import heapq
import math
from sympy import primerange

def smooth_via_heap(B: int, primes: [int]) -> [int]:
    """
    Generate all P-smooth numbers < B in ascending order
    using a priority queue.
    - B:     cutoff
    - primes: sorted list of primes [p1, p2, ..., pk]
    """
    heap = [1]
    seen = {1}
    result = []

    while heap:
        x = heapq.heappop(heap)
        if x > B:
            break
        result.append(x)
        for p in primes:
            next = x * p
            if next < B and next not in seen:
                seen.add(next)
                heapq.heappush(heap, next)
    return result

def valuation(N: int, P: int) -> int:
    """
    Compute the P-adic valuation of N.

    The P-adic valuation v_P(N) is the largest exponent k such that P**k divides N.
    """
    result = 0
    while (N%P) == 0:
        N //= P
        result += 1
    return result

def main(P: int, N: int, T: int, ilp: bool = False):
    primes = list(primerange(2, P+1))
    smooth = smooth_via_heap(P*T, primes) # worst-case we want to use a power of P that is >=T
    budget = {p: 0 for p in primes}
    target = {s: 0 for s in smooth}

    for i in range(1, N+1):
        t = i
        # Divide out the small primes
        for p in primes:
            while (t % p) == 0:
                t //= p
                budget[p] += 1

        # Figure out what factor is needed now
        for s in smooth:
            if t*s >= T:
                target[s] += 1
                break

    # No objective function, just feasibility
    print( "max: ;" )

    # Earth-moving constraints
    cdf = 0
    lhs = []
    for s in reversed(sorted(target.keys())):
        lhs.append(f"a{s}")
        if target[s] == 0:
            continue
        cdf += target[s]
        print("em", s, ": ", "+".join(lhs), ">=", cdf, ";", sep="")

    # Prime budgets
    for p in primes:
        lhs = []
        for s in smooth:
            v = valuation(s,p)
            if v == 0:
                continue
            elif v == 1:
                lhs.append( f"a{s}" )
            else:
                lhs.append( f"{v}a{s}" )
        if lhs:
            print("p", p, ": ", "+".join(lhs), "<=", budget[p], ";", sep="")

    # Require integer variables?
    if ilp:
        print("int ", ",".join([f"a{s}" for s in smooth]), ";", sep="")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
            prog='rearrange_lp',
            description='Generate an lp-format program for rearranging small prime factors.')
    parser.add_argument('P', type=int)
    parser.add_argument('N', type=int)
    parser.add_argument('T', type=int, nargs='?', default=0)
    parser.add_argument('--ilp', action='store_true', help='Solve optimization problem as integer program')

    args = parser.parse_args()

    P = args.P
    N = args.N
    T = args.T or math.ceil(N/4)

    main(P, N, T, args.ilp)
