import math

# A non-rigorous evaluation of c_0, by explicitly computing $\lfloor y \rfloor \log(\lceil y/e \rceil / (y/e)) dy/y^2$ on pieces up to some 
# threshold $b$, then using e/2b as an estimate for the remainder (the rationale being that for large $y$ one has $\log(\lceil y/e \rceil / (y/e)) \approx \{ y/e\}/(y/e)$, that $\lfloor y \rfloor \approx y$, and $latex \{y/e\}$ is approximately $1/2$ on the average, so that the tail is heuristically $\approx \int_b^\infty y \frac{1}{2} / (y/e) dy/y^2 = e/2b$)



def compute_c0():
    sum = math.pi**2 / 6

    for k in range(1,200000):
        sum += (math.log(1 + 1/k))**2 - 1/k**2
        if k % 1000 == 0:
            print(f"After {k} iterations, initial sum is known to be in: [{(sum-1/(2*k*k))/(2*math.e)},{sum/(2*math.e)}]")

    integral = sum/(2*math.e) + 2 * math.exp(-2) - math.log(2) / (2 * math.e)

    a = math.e
    N = 2
    M = 2

    for i in range(50000):
        if N+1 < math.e * M:
            b = N+1
            N_new = N+1
            M_new = M
        else:
            b = math.e * M
            N_new = N
            M_new = M+1
        integral += N * ((-math.log(M) + math.log(b))/b - (-math.log(M) + math.log(a))/a) / math.e
        integral -= ((math.log(M)+1)*(math.log(b)-math.log(a)) - (math.log(b)**2 - math.log(a)**2)/2) / math.e

        if (i%1000) == 0:
            print(f"Integral in: [{integral - 1/(2*b*b)}, {integral}], possibly {integral - 1/(8*b*b)}")
        N = N_new
        M = M_new
        a = b
    return integral

def integrand(x):
    if x == 0:
        return 0
    return math.floor(1/x) * math.log((math.e * x * math.ceil(1/(math.e * x)))) * math.log(1/x)

# Integral of \lfloor 1/x \rfloor log ( ex \lceil 1/ex \rceil ) log(1/x) dx from 0 to 1
def compute_c1p():
    N = 550000
    sum = 0
    altsum = 0
    for k in range(1, N+1):
        x = k / N
        sum += integrand(x) / N
        altsum += math.floor(1/x) * math.log((math.e * x * math.ceil(1/(math.e * x)))) / N
    print(f"c0 is approximately {altsum/math.e}")
    return sum / math.e

# \sum_{k=1}^\infty \frac{1}{k}  \log\left( \frac{e}{k} \left\lceil \frac{k}{e} \right\rceil \right) 
def compute_c1pp():
    sum = 0
    for k in range(1, 100000):
        sum += math.log( (math.e/k) * math.ceil(k/math.e) ) / k
        if k % 1000 == 0:
            print(f"After {k} iterations, c1'' has summed to {sum}")
    return sum 

c0 = compute_c0()
c1p = compute_c1p()
c1pp = compute_c1pp()

# c'_1 + c_0 c''_1 - ec_0^2/2
c1 = c1p + c0 * c1pp - (math.e * c0**2) / 2
print(f"c0 = {c0}, c1p = {c1p}, c1pp = {c1pp}, c1 = {c1}")

#compute_c0()