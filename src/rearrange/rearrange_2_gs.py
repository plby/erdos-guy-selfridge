def w(N):
    return bin(N).count("1") # Or in Python 3.8+: N.bit_count()

def check(N,T):
    return N >= 2 * ((8*T + 2 - 2*w(T)) // 3)

N = 1
T = 1
while True:
    while check(N,T):
        T += 1
    print(2,N,T)
    N += 1
