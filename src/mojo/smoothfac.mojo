
from collections import Dict, List
from gurobi import Environment, Model
from math import ceil, ceildiv, exp, floor, fma, log, sqrt
from memory import ArcPointer
from os.env import getenv
from sys import argv
from testing import assert_true
from time import perf_counter_ns
from utils.numerics import isnan, nan


@always_inline
fn sum[T: DType](read data: List[Scalar[T]]) -> Scalar[T]:
    var s: Scalar[T] = 0
    for x in data:
        s += x[]
    return s


fn coliter[func: fn(Int, Int) capturing [_] -> None](f: List[Int32], _j: Int) raises:
    """Iterate over prime divisors of factor j."""
    var j: Int = _j
    while j > 1:
        var i: Int = Int(f[j])
        assert_true(i > 1)
        var fij: Int = 0
        while f[j] == i:
            assert_true((j//i)*i == j)
            fij += 1
            j = j//Int(i)
        func(Int(i), fij)


fn colreduce[T: DType, //, func: fn(Scalar[T], Int, Int) capturing [_] -> Scalar[T]](f: List[Int32], _j: Int, _r: Scalar[T]) raises -> Scalar[T]:
    """Reduce prime divisors of factor j."""
    var j: Int = _j
    var r: Scalar[T] = _r
    while j > 1:
        var i: Int = Int(f[j])
        assert_true(i > 1)
        var fij: Int = 0
        while f[j] == i:
            assert_true((j//i)*i == j)
            fij += 1
            j = j//Int(i)
        r = func(r, Int(i), fij)
    return r


fn sieve_smooth(N: Int, sqN: Int, mut f: List[Int32], mut p: List[Int32], mut px: List[Int32], mut c: List[Int64]) raises:
    """Sieve small prime divisors <= √N."""
    print('Sieving small primes... ', end='')
    var t0 = perf_counter_ns()
    for i in range(2, sqN+1):
        if f[i] != 1:
            continue  # i is not prime
        ci = 0
        j = i
        while j <= N:
            if f[j] == 1:
                f[j] = i # store smallest prime factor of j
            fj = j
            while fj > 1:
                q = fj//i
                if q*i != fj:
                    break
                fj = q
                ci += 1
            j += i
        p.append(i)
        px[i] = len(p)-1
        c.append(ci)
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
    print('Number of small primes:', len(p))
    print('Count of small divisors:', sum(c))


fn sieve_nonsmooth(N: Int, sqN: Int, mut f: List[Int32]) raises -> Int:
    """Sieve large prime divisors > √N."""
    print('Sieving large primes... ', end='')
    var t0 = perf_counter_ns()
    var nf_nsmth = 0
    for i in range(sqN+1, N+1):
        if f[i] != 1:
            continue  # i is not prime
        ci = 0
        j = i
        while j <= N:
            f[j] = 0 # mark as non-smooth
            fj = j
            while fj > 1:
                q, m = divmod(fj, i)
                if m != 0:
                    break
                fj = q
                ci += 1
            j += i
        f[i] = -Int32(ci)  # store factor count (negate to differentiate from smooth factors)
        nf_nsmth += ci
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
    return nf_nsmth


fn deflate(N: Int, T: Int, sqN: Int, mut f: List[Int32], p: List[Int32], px: List[Int32], c: List[Int64], mut d: List[Float64]) raises:
    """Deflate small prime counts."""
    print('Deflating small prime counts... ', end='')
    var t0 = perf_counter_ns()
    # init d with c
    for ix in range(len(p)):
        d[ix] = Float64(c[ix])
    for j in range(sqN+1, N+1):
        if f[j] >= 0:
            continue  # j is not prime
        # clear 'greedy complements' (from previous iterations)
        l = 2*j
        while l <= N:
            assert_true(f[l] <= 0)
            f[l] = 0 # mark as non-smooth
            l += j
        # 'greedy complement' of j
        l = Int(ceildiv(T, j))
        assert_true(l <= sqN)
        assert_true((j<T) != (l==1))
        if l > 1:
            f[j*l] = -Int32(l)  # store 'greedy complement' (negate to differentiate from smooth factors)
        # 'deflate' divisors of greedy complement
        var cj = Float64(-f[j])
        while l > 1:
            k = Int(f[l])
            assert_true(k > 0 and (l//k)*k == l)
            d[px[k]] -= cj
            # print(j, l, k, d[k], cj, c[j])
            assert_true(d[px[k]] >= 0.0)
            l /= k
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')


fn lp_add_col(mod: ArcPointer[Model], f: List[Int32], px: List[Int32], j: Int, mut idx: List[Int32], mut val: List[Float64]) raises:
    idx.clear()
    val.clear()
    @parameter
    fn add_coef(i: Int, fij: Int):
        idx.append(px[i])
        val.append(Float64(fij))
    coliter[add_coef](f, j)
    mod[].add_var(nnz=len(idx), idx=idx, val=val, obj=1.0, lb=0.0, ub=Float64('inf'), vtype=ord('C'), name='x' + String(j))

fn lp_price(f: List[Int32], px: List[Int32], j: Int, a: List[Float64]) raises -> Float64:
    @parameter
    fn rc(aj: Float64, i: Int, fij: Int) -> Float64:
        return aj - a[px[i]]*Float64(fij)
    return colreduce[rc](f, j, 1.0)

fn solve_lp(N: Int, T: Int, sqN: Int, f: List[Int32], p: List[Int32], px: List[Int32], d: List[Float64], mut cols: List[Int], mut x: Dict[Int,Float64], mut a: List[Float64], lp_method: Int) raises -> Float64:
    # setup LP
    env = Environment.create()
    # env.set_output_flag(0)  # disable all output
    env[].set_threads(1)
    env[].set_method(lp_method)
    print('Setting up initial LP... ', end='')
    t0 = perf_counter_ns()
    # create model
    mod = Model.create(env, 'smoothfac')
    # objective
    mod[].set_model_sense(-1)  # min: 1, max: -1
    # scaling factor for RHS
    var dmin: Float64 = Float64('inf')
    for ix in range(len(p)):
        dmin = min(dmin, d[ix])
    var dscal: Float64 = Float64(1 << Int(floor(log(dmin)/log(2.0))))
    # print('\n*** RHS scaling:', dscal)
    # constraints
    for ix in range(len(p)):
        mod[].add_constr(sense=ord('<'), rhs=d[ix]/dscal, name='A'+String(p[ix]))
    mod[].update()
    # variables
    ncols = 0
    cols.clear()
    idx = List[Int32]()
    val = List[Float64]()
    # add initial batch of columns
    for j in range(T, N+1):
        if f[j] <= 0:
            j += 1
            continue  # skip non-smooth factors
        assert_true(f[j] > 1)
        # print('*** add col j =', j)
        lp_add_col(mod, f, px, j, idx, val)
        cols.append(j)
        ncols += 1
        if ncols >= max(2*len(p), 1000):
            break
    mod[].update()
    # mod.write("smoothfac.mps")
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')

    x.clear()
    var objval: Float64
    while True:
        # solve LP
        mod[].optimize()
        assert_true(mod[].get_status() == 2)  # 2: optimal
        objval = mod[].get_objval()
        # ***DEBUG***
        # print('ub =', mod[].get_objval())
        # for i in range(2, sqN+1):
        #     if c[i] == 0.0:
        #         continue
        #     print('i =', i, 'a =', mod[].get_pi(row_idx[i]), 'c =', c[i])
        # j = T
        # for ix in range(min(len(cols), 5*sqN+1)):
        #     jx = cols[ix]
        #     while j < jx:
        #         assert_true(f[j] <= 0)
        #         print(' ', j, ' ', -f[j])
        #         j += 1
        #     print(('*' if mod[].get_x(ix) > 0.0 else '·'), j, ' ', mod[].get_x(ix))
        #     j += 1
        # ***ENDEBUG***
        # get dual variables
        for ix in range(len(p)):
            a[ix] = mod[].get_pi(ix)
        # sift through next batch of columns
        print('Sifting columns... ', end='')
        t0 = perf_counter_ns()
        nbatch = 0
        for j in range(T, N+1):
            if f[j] <= 0:
                continue  # skip non-smooth factors
            # reduced cost of column j
            var aj: Float64 = lp_price(f, px, j, a)
            # TODO: better numeric estimate
            if aj/objval > 1e-10:
                # print('*** adding column j =', j, ' aj =', aj)
                if j in x:
                    print('*** WARNING: column', j, 'already in working set')
                else:
                    lp_add_col(mod, f, px, j, idx, val)
                    cols.append(j)
                    x[j] = 0.0
                    ncols += 1
                    nbatch += 1
            if nbatch >= 200:
                break
        print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
        if nbatch == 0:
            print('No more columns found.')
            break
        print('Reoptimizing with', nbatch, 'columns added')
        mod[].update()
    
    # recheck columns
    # for j in range(T, N+1):
    #     if f[j] <= 0:
    #         j += 1
    #         continue  # skip non-smooth factors
    #     var aj = lp_price(f, px, j, a)
    #     if aj > 1e-8:
    #         var is_col: Bool = False
    #         try:
    #             is_col = cols.index(j) >= 0
    #         except ValueError:
    #             pass
    #         print('*** ERROR: j =', j, ' aj =', aj, '*' if is_col else '')

    # get solution from LP model
    x.clear()
    for jx in range(ncols):
        if mod[].get_vbasis(jx) != 0:  # skip non-basic variables
            assert_true(mod[].get_x(jx) == 0.0)
            continue
        x[cols[jx]] = mod[].get_x(jx)*dscal


    return mod[].get_objval()*dscal


fn compute_ub(N: Int, T: Int, sqN: Int, f: List[Int32], p: List[Int32], px: List[Int32], c: List[Int64], a: List[Float64]) raises -> Float64:
    """Compute upper bound for full LP."""
    print('Computing upper bound... ', end='')
    var t0 = perf_counter_ns()
    var ub = 0.0
    # duals of small primes (from LP)
    for ix in range(len(p)):
        ub += a[ix]*Float64(c[ix])
    # duals of large primes
    for j in range(T, N+1):
        cj = -Int(f[j])
        if cj <= 0:
            continue
        i, m = divmod(j, cj) if cj > 1 else (j, 0)
        if m != 0:
            i = j
        if i < j:
            assert_true(f[i] < 0)
            cj = -Int(f[i])
        # print(j, i, cj, sqN, T)
        assert_true(i > sqN)
        aj = 0.0
        l = j
        while l <= N:
            # print('   ', l, l//i)
            var al: Float64 = lp_price(f, px, l//i, a)
            aj = max(aj, al)
            l += i
        ub += aj*cj
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
    return ub


fn greedy(N: Int, T: Int, sqN: Int, f: List[Int32], p: List[Int32], px: List[Int32], d: List[Float64], mut x: Dict[Int,Float64], mut g: List[Int]) raises -> Int:
    """Final greedy step using residual prime divisors of LP solution."""
    print('Computing greedy factors... ', end='')
    var t0 = perf_counter_ns()
    # 'deflate' LP factors
    for i in range(len(p)):
        g[i] = Int(d[i])
    for e in x.items():
        var j: Int = e[].key
        var xj: Int = Int(e[].value)
        if xj == 0:
            continue
        l = j
        while l > 1:
            i = Int(f[l])
            assert_true(i > 1 and (l//i)*i == l)
            g[px[i]] -= xj
            assert_true(g[px[i]] >= 0)
            l = l//i
    # greedy
    var nf_grdy: Int = 0
    var ix: Int = len(p)-1
    while ix >= 0:
        if g[ix] == 0:
            ix -= 1
            continue
        # 'greedy complement' of i
        m = 0
        var pi: Int = Int(p[ix])
        var j: Int = pi*ceildiv(T, pi)
        while j <= N:
            if f[j] <= 0:
                j += 1
                continue
            @parameter
            fn mindeg(m: Int64, i: Int, fij: Int) -> Int64:
                return min(m, g[px[i]]//fij)
            m = Int(colreduce[mindeg](f, j, Int.MAX))
            if m >= 1:
                break
            j += pi
        if j > N:
            break
        # record greedy factor
        nf_grdy += m
        x[j] = x.setdefault(j, 0.0) + 1.0
        # 'deflate' LP factors
        while j > 1:
            k = Int(f[j])
            assert_true(k > 1 and (j//k)*k == j)
            g[px[k]] -= m
            assert_true(g[px[k]] >= 0)
            j = j//k
    # assert sum(sc.values()) == nf_lp + nf_grdy
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
    return nf_grdy


fn main() raises:

    var args = argv()
    if len(args) < 2:
        print('Error: missing argment \'N\'')
        return

    # fractional mode (write LP solution into factorization file)
    var fractional = getenv('SMTHFC_FRACTIONAL') == '1'
    # LP solution method (0: primal simplex, 1: dual simplex, 2: interior point method)
    var lp_method = Int(getenv('SMTHFC_LP_METHOD', '-1'))
    # fixed threshold (disable OEIS mode)
    var T_fix = Int(getenv('SMTHFC_FIX_THRES', '-1'))

    var N: Int = Int(args[1])
    var sqN: Int = Int(ceil(sqrt(Float64(N))))

    print('Input: N =', N, ' √N =', sqN)

    # allocate memory for factorization
    print('Allocating ' + String(4*(N+1)//1024//1024) + 'MB memory... ', end='')
    t0 = perf_counter_ns()
    # f value encoding
    #  - f[j] > 1: j is smooth with smallest prime f[j]
    #  - f[j] = 1: j was not visited
    #  - f[j] = 0: j is not smooth
    #  - f[j] < 0: -f[j] is count of non-smooth factor j>=T (j=i*ceil(N/i) with i prime and ceil(N/i)<=√N)
    var f: List[Int32] = List[Int32](length=N+1, fill=1)          # factorizations of smooth factors / counts of non-smooth factors
    var c: List[Int64] = List[Int64](capacity=sqN+1)              # prime divisor counts of N! for small primes < √N
    var p: List[Int32] = List[Int32](capacity=sqN+1)              # prime index of small primes (p[i] is i-th prime)
    var px: List[Int32] = List[Int32](length=sqN+1, fill=-1)      # inverse prime map (px[i] is index of prime i, or -1 if i is not prime)
    # additional LP storage
    var d: List[Float64] = List[Float64](length=sqN+1, fill=0.0)  # deflated prime divisor counts of N! for small primes < √N (RHS of LP)
    var a: List[Float64] = List[Float64](length=sqN+1, fill=0.0)  # optimal dual multipliers for small primes < √N
    var cols: List[Int] = List[Int](capacity=sqN+1)               # column mapping (LP column index -> factor j)
    var x: Dict[Int,Float64] = Dict[Int,Float64]()                # optimal values of LP (plus greedy values)
    # additional greedy workspace
    var g: List[Int] = List[Int](length=sqN+1, fill=0)            # deflated prime divisor counts of N! for small primes after LP solution
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')

    # sieve
    sieve_smooth(N, sqN, f, p, px, c)
    var nf_nsmth: Int = sieve_nonsmooth(N, sqN, f)
    print('Number of non-smooth factors: ' + String(nf_nsmth))

#    if N > 300_000:
#        return

    # print('--------------')
    # for j in range(2, N+1):
    #     print(j, f[j])

    var T_lo_ini: Int
    var T_hi_ini: Int
    if T_fix >= 0:
        T_lo_ini = 0
        T_hi_ini = Int.MAX
    elif N >= 300_000:
        # lower bound from eq. (1.8), upper bound heuristic
        T_lo_ini = Int(Float64(N)*(1/exp(1.0) - 0.3044190/log(Float64(N)) - 0.7555/(log(Float64(N))**2)))
        T_hi_ini = Int(ceil(1.01*T_lo_ini))
    elif N >= 50_000:  # some safety of true value 43632
        T_lo_ini = Int(ceildiv(Float64(N), 3.0))
        T_hi_ini = 2*T_lo_ini
    else:
        T_lo_ini = ceildiv(2*N, 7)
        if N == 56:
            T_lo_ini -= 1
        T_hi_ini = 2*T_lo_ini


    # bisect for lower bound
    var T_lo = T_lo_ini
    var T_hi = T_hi_ini
    var T_ub_lo: Int = 0
    var T_ub_hi: Int = Int.MAX
    while T_hi-T_lo > 1:

        var T: Int = (T_lo+T_hi)//2 if T_fix == -1 else T_fix
        print('Solving problem: N =', N, ' T =', T, ' √N =', sqN)
        print('*** N =', N, ' T_lo = ', T_lo, ' T_hi =', T_hi, ' T =', T)

        # deflate small prime counts
        deflate(N, T, sqN, f, p, px, c, d)

        # solve LP
        lp_val = solve_lp(N, T, sqN, f, p, px, d, cols, x, a, lp_method)

        objval = lp_val + nf_nsmth
        print('LP optimal value: ' + String(objval))

        if T > T_ub_lo:
            # compute upper bound
            var ub = compute_ub(N, T, sqN, f, p, px, c, a)
            print('Upper bound: ' + String(ub) + ' (gap: ' + String(round((ub-objval)/ub*10000.0)/100.0) + '%)')
            if ub < Float64(N):
                assert_true(T < T_ub_hi)
                T_ub_hi = T
            else:
                T_ub_lo = T

        # truncate and count LP factors
        var nf_lp: Int = 0
        for e in x.items():
            var j: Int = e[].key
            var xj: Int = Int(e[].value)  # TODO: beware of rounding errors
            if xj == 0:
                continue
            if not fractional:
                e[].value = Float64(xj)
            nf_lp += xj
        # assert_true(sum(sc.values()) == nf_lp)
        print('Number of LP factors: ' + String(nf_lp))

        var nf_grdy: Int = 0
        if not fractional:
            nf_grdy = greedy(N, T, sqN, f, p, px, d, x, g)
            print('Number of greedy factors:', nf_grdy)
        
        var nf_tot = nf_nsmth + nf_lp + nf_grdy
        print('Total factors:', nf_tot)

        if T_fix >= 0:
            break

        if nf_tot >= N:
            T_lo = T
        else:
            T_hi = T
    
    assert_true(T_fix >= 0 or T_lo > T_lo_ini)
    assert_true(T_fix >= 0 or T_hi < T_hi_ini)

    # bisect for upper bound
    assert_true(T_fix >= 0 or T_ub_lo > 0)
    assert_true(T_fix >= 0 or T_ub_hi < Int.MAX)
    while T_fix == -1 and T_ub_hi-T_ub_lo > 1:

        var T: Int = (T_ub_lo+T_ub_hi)//2
        print('Solving problem: N =', N, ' T =', T, ' √N =', sqN)
        print('*** N =', N, ' T_lo = ', T_ub_lo, ' T_hi =', T_ub_hi, ' T =', T)

        # deflate small prime counts
        deflate(N, T, sqN, f, p, px, c, d)

        # solve LP
        lp_val = solve_lp(N, T, sqN, f, p, px, d, cols, x, a, lp_method)

        objval = lp_val + nf_nsmth
        print('LP optimal value: ' + String(objval))

        # compute upper bound
        var ub = compute_ub(N, T, sqN, f, p, px, c, a)
        print('Upper bound: ' + String(ub) + ' (gap: ' + String(round((ub-objval)/ub*10000.0)/100.0) + '%)')

        if ub < Float64(N):
            T_ub_hi = T
        else:
            T_ub_lo = T

    if T_fix == -1:
        print('OEIS result: N =', N, ' t(N) =', T_lo, ' T_ub =', T_ub_hi)

    # write output file
    if len(args) > 2:
        print('Writing factorization file... ', end='')
        t0 = perf_counter_ns()
        with open(args[2], 'w') as file:
            file.write(N, ' ', T_lo, '\n')
            file.write("FACTORS\n")
            for j in range(T_lo, N+1):
                var xj: Float64 = x.get(j, nan[DType.float64]())
                if fractional and not isnan(xj):
                    file.write(j, ' ', xj, '\n')
                    continue
                var fj: Int
                if not isnan(xj):
                    fj = Int(xj)
                else:
                    fj = -Int(f[j])
                    if fj <= 0:
                        continue
                    i, m = divmod(j, fj) if fj > 1 else (j, 1)
                    if m == 0:
                        fj = -Int(f[i])
                file.write(j, ' ', fj, '\n')
            # residual = prod(pow(i, int(ci)) for (i, ci) in enumerate(c) if ci > 0.0)
            # if residual > 1:
            #     print(f"{residual} 1", file=file)
            print("CERTIFICATE", file=file)
            # TODO
        print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')

# vim: ts=4 sw=4 sts=4 et

