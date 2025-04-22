import math

data = [
    [100000, 33462, 33646],
    [200000, 67703, 67704],
    [300000, 101903, 101907],
    [400000, 136143, 136149],
    [500000, 170456, 170464],
    [600000, 204811, 204821],
    [700000, 239187, 239210],
    [800000, 273604, 273627],
    [900000, 308029, 308042],
    [1000000, 342505, 342508],
    [2000000, 687796, 687802],
    [3000000, 1033949, 1033966],
    [4000000, 1380625, 1380654],
    [5000000, 1727605, 1727609],
    [6000000, 2074962, 2074994],
    [7000000, 2422486, 2422540],
    [8000000, 2770212, 2771082],
    [9000000, 3118129, 3118161],
    [10000000, 3466235, 3466252],
    [20000000, 6952243, 6952266],
    [30000000, 10444441, 10444466],
    [40000000, 13940484, 13942660],
    [50000000, 17439282, 17439315],
    [60000000, 20940210, 20940249],
    [70000000, 24442818, 24442966],
    [80000000, 27946958, 27947028],
    [90000000, 31452431, 31452469],
    [100000000, 34958725, 34959067],
    [200000000, 70064782, 70065125],
    [300000000, 105218403, 105218496],
    [400000000, 140401212, 140401354],
    [500000000, 175605266, 175605553],
    [600000000, 210825848, 210825988],
    [700000000, 246059851, 246060084]
]

        
def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def logfac(N):
    sum = 0
    for i in range(1, N+1):
        sum += math.log(i)
    return sum

# Test if  \sum_{p > \frac{t}{\sqrt{t}+1}} \left\lfloor \frac{N}{p} \right\rfloor \log \left( \frac{p}{t} \left\lceil \frac{t}{p} \right\rceil \right) > \log N! - N \log t

def criterion( N, t ):
    sum = 0
    threshold = logfac(N) - N * math.log(t) + 0.0001 
    tfail = t/(math.floor(math.sqrt(t)))
    for p in reversed(range(1, N + 1)):
        if p < tfail:
            return False
        if is_prime(p):
            sum += math.floor(N/p) * math.log((p/t)*math.ceil(t/p))
            if sum > threshold:
                return True
    return False

def best_t( N, init ):
    t = init
    while not criterion(N, t):
        t += 10
    while criterion(N,t):
        t -= 1
    print(f"Testing {N}: needed to increment by {t-init}") 
    return t



c0 = 0.30441901
c1 = 0.7554808

def round_to_sigfigs(n, sigfigs):
    if n == 0:
        return 0
    import math
    digits = int(math.floor(math.log10(abs(n)))) + 1
    factor = 10**(max(digits - sigfigs,0))
    return (n // factor) * factor


for x in data:
    N, t_lower, t_upper = x
    diff = t_upper - t_lower

    print(N, t_lower, best_t(N,t_lower))
#    approx = int(N / math.e - c0 * N / math.log(N) - c1 * N / (math.log(N) ** 2))
#    print(f"$\\num{{{N}}}$ & $[\\num{{{t_lower}}}, \\num{{{t_upper}}}]$ & $\\num{{{diff}}}$ & $\\num{{{round_to_sigfigs(approx,7)}}}$ \\\\")
