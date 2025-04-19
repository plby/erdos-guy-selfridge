import math

# A non-rigorous evaluation of c_0, by explicitly computing $\lfloor y \rfloor \log(\lceil y/e \rceil / (y/e)) dy/y^2$ on pieces up to some 
# threshold $b$, then using e/2b as an estimate for the remainder (the rationale being that for large $y$ one has $\log(\lceil y/e \rceil / (y/e)) \approx \{ y/e\}/(y/e)$, that $\lfloor y \rfloor \approx y$, and $latex \{y/e\}$ is approximately $1/2$ on the average, so that the tail is heuristically $\approx \int_b^\infty y \frac{1}{2} / (y/e) dy/y^2 = e/2b$)


def integrand(x):
    return math.floor(x) * math.log(math.ceil(x/math.e) / (x/math.e)) / x**2

sum = math.pi**2 / 6

for k in range(1,200000):
    sum += (math.log(1 + 1/k))**2 - 1/k**2
    if k % 1000 == 0:
        print(f"Sum after {k} iterations: {sum/(2*math.e)}")


integral = sum/(2*math.e) + 2 * math.exp(-2) - math.log(2) / (2 * math.e)
base = integral

print(integral)

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

    # We can note that the function we are integrating is lower bounded by 0 and upper bounded by e/2 *(1 + 1/(1 + ex))
#    print(f"After incorporating [{a},{b}], the estimate for $c_0$ is {(integral + math.exp(1)/(2*b))/math.e}.") 
#    print(f"Bounds: [{integral/math.e}, {(integral + (math.exp(1)/b + math.log(1 + math.exp(1)/b))/2)/math.e}]")
#    print(f"Running upper bound: {integral2}")
    print(f"Integral in: [{integral - 1/(2*b*b)}, {integral}], possibly {integral - 1/(8*b*b)}")
    N = N_new
    M = M_new
    a = b

print(base-integral)

# $$ 0.304419004 \leq c_0 \leq 0.304419017.$$

