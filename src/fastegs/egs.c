#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <memory.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <primecount.h>         // uncomment if you want to use primecount_pi rather than loading pi.dat
#include <primesieve.h>
#include <primesieve/iterator.h>

#define MAXPI   ((1<<16)-1)     // cannot exceed 2^16
#define MAXP    821603          // MAXPI-th prime
#define MAXM    ((1<<24)-1)     // must exceed sqrt(t), ideally at least t^(5/8)
#define MAXN    ((1L<<40)-1)    // limited by pi-value range in pi.h, which could easily be increased
#define MAXF    59804712        // sum_m (1+#{p|m}) over MAXP-smooth m <= MAXM

/*
    The tables P,PI,F,M are independent of N <= 2^40 and computed at startup (about 0.5s)
    and are read-only after that.  They use about 300MB memory (this is most of what we will use)
*/

static int32_t *P;              // P[n] is the nth prime for n up to MAXPI, and we put P[0] = 1
static int32_t *PI;             // PI[n] = pi(n) for n <= MAXP (in particular, P[PI[p]]=p for primes p)

// Factorizations of cofactors m are zero-terminated lists of pp's
// We only consider m that are MAXP-smooth, so pi fits in 16-bits (could extend to 24 bits and make e 8bits)
static struct pp {
    uint16_t pi;                // index into P
    uint16_t e;
} *F;                           // concatenation of zero-terminated factorizations in descending order by pi
int32_t *M;                     // F[M[m]] holds the factorization of m <= MAXM for MAXP-smooth m, M[m]=0 ow


static inline double get_time (void) // accurate to at least 10ms
    { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (double)t.tv_sec + (double)t.tv_nsec / 1000000000.0L; }

static inline int64_t cdiv (int64_t a, int64_t b)     // ceil(a/b)
    { return (a+b-1) / b; }

static inline int v2 (uint32_t x) { return __builtin_ctz(x); }

static inline int64_t pi (int64_t n)
    { return n <= MAXP ? PI[n] : primecount_pi(n); }

static inline int64_t max (int64_t a, int64_t b)
    { return a<b?b:a; }

static inline int64_t min (int64_t a, int64_t b)
    { return a>b?b:a; }

static inline primesieve_iterator primesieve_start(int64_t minp, int64_t maxp)
{
    primesieve_iterator ctx;
    primesieve_init (&ctx);
    primesieve_jump_to (&ctx, minp, maxp);
    return ctx;
}

static inline void primesieve_stop (primesieve_iterator *ctx)
    { primesieve_free_iterator (ctx); }

static inline int is_prime (int64_t n)
    { return pi(n) > pi(n-1); } // TODO: replace this with BPSW primality test

// setup here is independent of N and t and will handle any N <= MAXN = 2^40
// this takes about 0.5s and uses around 200MB which is neglible if we are
// going to be running computations for many t and N
void general_setup (void)
{
    // compute the list of primes p <= MAXP, storing the nth prime in P[n] (we set P[0]=1 for convenience)
    P = malloc((MAXPI+1)*sizeof(*P)); assert(P);  P[0] = 1;
    PI = calloc((MAXP+1),sizeof(*PI)); assert (PI);
    primesieve_iterator ctx = primesieve_start (0,MAXP);
    for ( int32_t p = 0, n = 1 ; (p = primesieve_next_prime(&ctx)) <= MAXP ; ) { P[n] = p; PI[p] = n++; }
    primesieve_stop(&ctx);
    assert (P[MAXPI] == MAXP);

    // set Pi[n] to pi(n) for all n (using values already set for n prime)    
    for ( int32_t n = 1 ; n <= MAXP ; n++ ) if ( !PI[n] ) PI[n] = PI[n-1];

    // set M[m]=n where p_n is the largest odd prime divisor of M up to MAXP (and zero if no such p_n exists)
    M = calloc(MAXM+1,sizeof(*M));
    for ( int32_t pi = 2 ; pi <= MAXPI ; pi++ ) for ( int p=P[pi], q = p ; q <= MAXM ; q+=p ) M[q] = pi;

    // Compute F and update M so that F[M[m]] holds the factorization of m for all MAXP-smooth m <= MAXM
    // Factorizations are zero terminated lists of pp in reverse order by prime
    F = malloc(2*(MAXF+1)*sizeof(*F));
    struct pp *f = F; (f++)->pi = 0;                // skip offset 0
    for ( int32_t m = MAXM ; m > 1 ; m-= 2 ) {      // handle odd m first
        struct pp *g = f;
        g->pi = M[m]; g->e = 0;
        int32_t q = m;
        for ( ; M[q] ; q /= P[M[q]] ) {
            if ( M[q] == g->pi ) { g->e++; }
            else { (++g)->pi = M[q]; g->e = 1; }
        }
        if (q!=1) { M[m] = 0; continue; }
        g++; (g++)->pi = 0;
        M[m] = f-F; f = g;
    }
    M[1] = f-F; (f++)->pi = 0;
    for ( int32_t m = MAXM-1 ; m > 1 ; m-= 2 ) {    // now handle even m
        int32_t e = v2(m), q = m>>e;
        if ( !M[q] ) { M[m] = 0; continue; }
        struct pp *g = F+M[q];
        M[m] = f-F;
        while (g->pi) *f++ = *g++;
        f->pi = 1; (f++)->e = e; (f++)->pi = 0;
    }
    assert (f-F == MAXF+1);
}

static inline int64_t fcnt (int64_t *E, int64_t e, struct pp *f) 
{   // computes min(e,v_m(P^E)) where m has factorization f and v_m(n) is largest r for which m^r|n
    for ( ; f->pi ; f++ ) e = min(e,E[f->pi]/f->e);
    return e;
}

// tfac returns a lower bound on the number of factors in any factorization of N! into parts of size >= t
// we require N/3 < t < N/2 (this simplifies the algorithm and will always hold in the context of checking EGS
// When check_feasibility is nonzero the algorithm will return 0 if discovers that the greedy algorithm cannot
// produce >= N factors (this check is made before any suboptimal cofactors are used), and when this occurs
// it applies not only to the inputs (N,t) but also to inputs (N,t') for all t' > t
int64_t tfac (int64_t N, int64_t t, int feasible, int verbosity)
{
    double start = get_time();
    if ( verbosity > 2 ) printf("tfac(%ld,%ld)\n",N,t);
    assert (N >= 10 && N < MAXN && 3*t >= N && 2*t < N);
    int32_t sqrtN = (int32_t)sqrt(N);
    int32_t s = sqrt(t); assert(s*(s-1) < t); while ( (int64_t)s*(s-1) < t ) s++;
    assert (s <= MAXP);
    int32_t maxpi = PI[s-1];
    
    // Compute p-adic valuations of N! for p < s by setting E[pi(p)] = v_p(N!)
    // We could cache and resue this but it takes negligible time to recompute it (under 1ms for N <= 10^12)
    int64_t *E = calloc(maxpi+1,sizeof(*E));
    for ( int32_t i = 1 ; i <= maxpi ; i++ ) for ( int64_t q = P[i] ; q <= N ; q *= P[i] ) E[i] += N/q;

    // Compute upper bound maxm on the smooth cofactors m we will use.
    // We want the optimal m = cdiv(t,p) for p to be (p-1)-smooth (this simplifies the algorithm)
    // maxm should be t^delta for some 1/2 < delta < 2/3
    // Larger values of delta use more memory/time but give slightly better factorizations
    // heuristically delta = 5/8 seems close to optimal for speed
    int32_t maxm = min((int32_t)pow(cdiv(N,3),0.625),MAXM);
    int32_t *Ms = malloc((maxm+1)*sizeof(*Ms)); Ms[0] = 0;
    for ( int32_t m = 1, *p = Ms+1 ; m < s ; m++ ) *p++ = m;
    int32_t numm = s;
    for ( int32_t m = s ; m <= maxm ; m++ ) { if ( F[M[m]].pi && F[M[m]].pi <= PI[t/m] ) Ms[numm++] = m; }
    Ms = realloc(Ms, numm*sizeof(*Ms));
    numm--; maxm = Ms[numm];

    /*
      We begin by constructing factors m*p >= t with m minimal for all p >= s >= sqrt(t).
      In this range we can always make m = cdiv(t,p) optimal for p.
      For each prime p in [s,N] we will have n=v_p(N!) identical factors m*p.
      We don't really care about p, we just need the values of m and n.
      We will treat m = ceil(t/s),...,1 in descending order.  When m is large, p will
      be small, and the value of n will be changing rapidly, so we just enumerate primes p
      and compute update m and n as we go.  But when m is small it is much more
      efficient to determine the exact range of p applicable to each (m,n) pair and
      then count the primes in this range, rather than enumerate them (we can us a
      precomputed table of pi(x) values to do this very quickly).

      The crossover point between "large" and "small" m is a performance parameter that
      depends on the relative speed of primesieve vs primecount.
      t^(1/5) seems like a reasonable heuristic value for large t
    */

    int64_t m = cdiv(t,s);                                  // largest possible m for p >= s
    assert (m <= maxm && Ms[m] == m);
    int64_t mid = min((int64_t)pow(t,0.2),(t-1)/sqrtN);     // m > mid are large, m <= mid are small
    if ( (int64_t)sqrtN*mid >= t ) mid = (t-1)/sqrtN;       // force p > sqrtN for m < mid (handy)

    if ( verbosity > 2 ) printf ("N=%ld, t=%ld, sqrt(N)=%d, s=%d, maxpi=%d, maxm=%d, numm=%d, mid=%ld (%.6fs)\n", N, t, sqrtN, s, maxpi, maxm, numm, mid, get_time()-start);

    primesieve_iterator ctx = primesieve_start(s,(t-1)/mid);
    int64_t p; // p will fit in 32 bits but we want to mults at 64-bits
    int64_t cnt = 0; // running count of the number of factors (goal is to get this >= N)

    // handle primes in [s,sqrt(N)] (here we are happy to recompute n for each p, this phase takes no time)
    while ( (p = primesieve_next_prime(&ctx)) <= sqrtN) {
        while ( (m-1)*p >= t ) m--;     // minimize m with m*p >= t
        int64_t n = N / p + N / (p*p);  // compute n = v_p(N!) (expensive but negligible)
        // update valuations in E to reflect our use of n cofactors m
        for ( struct pp *f = F+M[m] ; f->pi ; f++ ) E[f->pi] -= n*f->e;
        cnt += n;
    }

    if ( verbosity > 2 ) printf("cnt=%ld for p in [s,sqrt(N)], m=%ld (%.6fs)\n", cnt, m, get_time()-start);

    int64_t pmmax = (t-1) / (m-1);  // largest p for this m
    assert (p>pmmax || m==cdiv(t,p));
    int64_t n = N / (sqrtN+1);      // smallest n for p > sqrtN
    int64_t pnmax = N / n;          // largest p for this n
    int64_t plmmax = (t-1) / mid;   // largest p for m > mid

    // handle primes in (sqrtN,plmmax] with large m using primesieve (this should take about half the time if mid is optimal)
    for ( ; p <= plmmax ;) {
        while ( p > pmmax ) { m--; pmmax = (t-1)/(m-1); }   // update m
        while ( p > pnmax ) { n--; pnmax = N/n; }           // update n
        int64_t pmax = min(pmmax,pnmax); int64_t x = 0;
        for ( x += n ; (p = primesieve_next_prime(&ctx)) <= pmax ; x+=n );
        for ( struct pp *q = F+M[m] ; q->pi ; q++ ) E[q->pi] -= x*q->e;
        cnt += x;
    }
    primesieve_stop(&ctx);
    int64_t lastpi = pi(plmmax);
    if ( verbosity > 2 ) printf("cnt=%ld for %ld p >= s with m < mid (%.6fs)\n", cnt, lastpi-maxpi, get_time()-start);

    // handle primes in (plmax,t] with small m in [mid,2] using primecount (this should take about half the time if mid is optimal)
    // here we iterate over m rather than p
    for ( m = mid ; m > 1 ; m-- ) {
        int64_t pmin = cdiv(t,m), pmax = (t-1)/(m-1);       // p in [pmin,pmax] are the p for this m
        n = N/pmin; pnmax = min(N/n,pmax);                  // p in [pmin,pnmax] are the p for this n
        while ( pmin <= pmax ) {
            int64_t nextpi = pi(pnmax);
            int64_t c = n*(nextpi-lastpi); cnt += c;        // number of factors for this m and n
            for ( struct pp *q = F+M[m] ; q->pi ; q++ ) E[q->pi] -= c*q->e;
            pmin = pnmax+1; n--; pnmax = min(N/n,pmax);     // proceed to the next n
            lastpi = nextpi;
        }
    }
    if ( verbosity > 2 ) printf("cnt=%ld for %ld p in [s,t) (%.6fs)\n", cnt, lastpi-maxpi, get_time()-start);

    // Finally, handle primes p  in [t,N].  Here m=1 and n=3,2,1 (3 only if N=3t and t is prime)
    int64_t nextpi = pi(N/2);
    cnt += (N == 3*t && is_prime(t) ? 1 : 0) + 2*(nextpi-lastpi);
    lastpi = nextpi; nextpi = pi(N);
    cnt += nextpi-lastpi;
    if ( verbosity > 2 ) printf("cnt=%ld for %ld p in [s,N] (%.6fs)\n", cnt, nextpi-maxpi, get_time()-start);

    // Verify our assumption that we can use the optimal m for all p >= s (this will be verified again at the end but
    // the check is cheap so we do it now before a possible feasibility check).
    for ( int32_t i = 1 ; i <= maxpi ; i++ ) assert(E[i] >= 0);

    /*
        At this point we have dealt with all prime factors p >= s of N! using the optimal m = cdiv(t,p)
        for each p and updated E so that

            P^E := prod_{1<=i<=pimax} P[i]^E[i]

        is the divisor of N! that we still need to factor.

        We now process primes p=P[i] in descending order with cofactors m=M[j] in ascending order.
        We need p*m >= t but no longer insist that m = cdiv(t,p) is minimal (because eventually
        this will not be possible), and to simplify the algorithm we require m to be (p-1)-smooth.
        This ensures that the sets of exponents in E associated to p and m are disjoint, and if
        m >= cdiv(t,p) is (p-1)-smooth, so is every m' > m, which simplifies the code.

        This deviates from the original greeyd approach but for large N we are happy to obtain a
        lower bound on t that is potentially slightly worse in return for a faster/simpler algorithm.
        For N >= 10^6 there is plenty of room to spare. Our goal is to efficiently handle N in [10^6,10^12].

        This phase takes about 1-2 percent of the time.
    */
    int32_t pimin = pi(cdiv(t,maxm))+1;
    for ( int32_t i = maxpi, j = cdiv(t,s) ; i >= pimin ; i-- ) {
        // Update j so that all m' >= m=sM[j] are valid for use with p=P[i]
        while ( (int64_t)P[i]*Ms[j] < t || F[M[Ms[j]]].pi >= i ) j++;
        if ( feasible && Ms[j] >= 2*P[i] ) break;
        struct pp *f = F+M[Ms[j]]; assert(f);
        int64_t e = fcnt(E,E[i],f);    
        if ( e < E[i] ) { // if we cannot completely remove p from P^e using m, try p^2 and m=cdiv(t,p^2)
            struct pp *g = F+M[cdiv(t,(int64_t)P[i]*P[i])];
            e = fcnt(E,E[i]/2,g);
            if ( e ) { cnt += e; for ( E[i] -= 2*e ; g->pi ; g++) E[g->pi] -= e*g->e; }
            e = fcnt(E,E[i],f); // recompute e (may have changed)
        }
        if ( e ) { cnt += e; for ( E[i] -= e ; f->pi ; f++ ) E[f->pi] -= e*f->e; }
        if ( E[i] ) { // p still divides P^E?
            // try a larger m
            e = 0;
            for ( int32_t k = j+1 ; k <= numm ; k++ ) {
                struct pp *g = F+M[Ms[k]];
                int64_t x = fcnt(E,E[i],g);
                if ( x > e ) {  e = x; f=g; if ( e == E[i] ) break; }
            }
            if ( e ) { cnt += e; for ( E[i] -= e ; f->pi ; f++ ) E[f->pi] -= e*f->e; }
            if ( E[i] ) { // p still divides P^E?
                // try a larger m for p^2, here we assume m=cdiv(t,p^2) < s so that Ms[m] = m
                e = 0; f = 0; m = cdiv(t,(int64_t)P[i]*P[i])+1; assert (Ms[m]==m);
                for ( int32_t k = m ; k <= numm ; k++ ) {
                    struct pp *g = F+M[Ms[k]];
                    int64_t x = fcnt(E,E[i]/2,g);
                    if ( x > e ) {  e = x; f=g; if ( e == E[i] ) break; }
                }
                if ( e ) { cnt += e; for ( E[i] -= 2*e ; f->pi ; f++ ) E[f->pi] -= e*f->e; }
                // we may still have E[i] > 0 but we will usually have E[i] <= 1
            }
        }
    }
    if ( feasible ) {
        long double ebits = 0, epsilon = 0.0000000000000001L; // 10^{-15} < 2^52
        for ( int32_t i = 1 ; i <= maxpi ; i++ ) ebits += E[i]*log(P[i]+epsilon); // make sure we get an upper bound
        return cnt + floorl(ebits/log(t-epsilon));
    }

    if ( verbosity > 2 ) printf("cnt=%ld after initial pass of p in (cdiv(t,maxm),s) (%.6fs)\n", cnt, get_time()-start);
    while ( maxpi && !E[maxpi] ) maxpi--;
    free(Ms);

    // Now we want to use up whatever is lest the best we can.  This should consist almost entirely of primes < t^1/3 and takes no time
    int64_t good = 5*cdiv(t,4);  // we will settle for factors in [t,good]
    struct pp c[16]; c->pi = 0; c->e = 0; // factorization of product of primes we are assembling (which may exceed maxm)
    while (maxpi) {
        while ( maxpi && !E[maxpi] ) maxpi--;
        if ( !maxpi ) break;
        int32_t i = maxpi;
        int64_t q = P[i];
        struct pp *f=c; f->pi = i; f->e = 1; (++f)->pi = 0;
        E[i]--; // update E as we go, we will undo the update if we get stuck
        while ( i && !E[i] ) i--;
        if ( !i ) break;
        while ( i && q*P[i] < good ) {
            q *= P[i]; E[i]--;
            if ( (f-1)->pi == i ) { (f-1)->e++; } else { f->pi = i; f->e = 1; (++f)->pi = 0; }
            while ( i && !E[i] ) i--;
        }
        if ( !i && q < t ) break;
        int64_t e = 1 + fcnt(E,E[c->pi]/c->e,c+1);    // 1+v_q(P^E) (we aleady removed one factor of q)
        if ( q < t ) {
            assert (q>s);
            int64_t b = 0; struct pp *g = 0;
            // first look for a cofactor smaller than the smallest prime divisor of q
            for ( m = cdiv(t,q) ; m < P[(f-1)->pi] ; m++ ) {
                int64_t x = fcnt(E,e,F+M[m]);
                if ( x > b ) { b = x; g = F+M[m]; }
                if ( x == e ) break;
            }
            if ( b ) {
                while ( g->pi ) { E[g->pi] -= g->e; *f++ = *g++; }
                f->pi = 0;
            } else {
                if (!i) break;
                q *= P[i]; E[i]--;  assert (q >= t);
                if ( (f-1)->pi == i ) { (f-1)->e++; } else { f->pi = i; f->e = 1; (++f)->pi = 0; }
                b = 1 + fcnt(E,E[c->pi]/c->e,c+1);
                assert(b);
            }
            e = b;
        }
        cnt += e--;
        for ( f = c ; f->pi ; f++ ) E[f->pi] -= e*f->e;
        c->pi = 0;
        maxpi = i;
    }
    for ( struct pp *f = c ; f->pi ; f++ ) E[f->pi] += f->e; // restore any patial factor we did not remove so we can report/check remainder
    while ( maxpi && !E[maxpi] ) maxpi--;
    int64_t q = 1; for ( int32_t i = 1 ; i <= maxpi ; i++ ) { assert(E[i] >= 0); for ( int32_t e = 0 ; e < E[i] ; e++ ) { q *= P[i]; assert (q<t); } }
    if ( verbosity > 2 ) printf("cnt=%ld after final pass, remainder is %ld (%.6fs)\n", cnt, q, get_time()-start);
    free (E); E = 0;
    return cnt;
}

// returns a value of t >= N/3 that yields a good lower bound on t(N), or 0 if no such t can be found
int64_t tbound (int64_t N, int optimal, int verbosity)
{
    int64_t t = cdiv(N,3);
    int64_t cnt = tfac(N,t,verbosity,0);
    if ( cnt < N ) return 0;
    int64_t tmin = t, tmax = (2*N)/5;

    /*
        We use a modified bisection search of [tmin,tmax) for t with tfac(N,t) >= N but tfac(N,t+1) < N that
        uses the excess/deficit (tfac(N,t)-N) to choose the next bisection point.
    */
    while ( tmin < tmax-1 ) {
        if ( cnt >= N ) tmin = max(t,tmin); else tmax = min(t,tmax);
        if ( verbosity > 1 ) fprintf (stderr,"t=%ld gave %ld extra factors, new t-range is [%ld,%ld)\n", t, cnt-N, tmin, tmax);
        t = round(exp(log(t)+(cnt-N)*log(t)/N));
        if ( t <= tmin ) t = max((3*tmin+tmax)/4,tmin+1);
        if ( t >= tmax ) t = min((tmin+3*tmax)/4,tmax-1);
        cnt = tfac(N,t,0,verbosity);
    }
    assert (tmax < (2*N)/5);
    if ( ! optimal ) return tmin;
    if ( verbosity > 0 ) printf("t(%ld) >= %ld proved\n", N,tmin);

    /* Now use a binary search to get an upper bound on the best possiible t that tfac(N,t) could return */
    int64_t low = tmin;
    int64_t high = (2*N)/5;
    cnt = tfac(N,high,1,verbosity);
    assert (cnt < N);
    while ( low < high-1 ) {
        int64_t mid = (low+high)/2;
        cnt = tfac(N,mid,1,verbosity);
        if ( cnt < N ) { high = mid; tmax = mid; } else { low = mid; }
    }
    assert (tmax > tmin);
    if ( verbosity > 0 ) printf("t(%ld) >= %ld cannot be proved\n",N,tmax);
    if ( verbosity > 0 ) printf("checking %ld values of t\n", tmax-tmin-1);
    for ( t = tmin+1 ; t < tmax ; t++ ) {
        if ( tfac(N,t,0,verbosity) > N ) {
                tmin = max(tmin,t);
                if ( verbosity >= 0 ) printf("\rt(%ld) >= %ld proved\r", N,tmin);
        }
    }
    return tmin;
}

static void usage (void)
{
    fprintf(stderr,
        "Usage: egs [-v level] [-h filename] [-c] [-o] N-range [t]\n"
        "       -v level      integer verbosity level -1 to 3 (optional, default is 0)\n"
        "       -h filename   hint-file with records N:t (required if range of N is specified)\n"
        "       -c            create hint-file rather than reading it (must be specified in combination with -h)\n"
        "       -o            use the best t for which the algorithm can prove t(N) >= t (optional)\n"
        "       N-range       integer N or range of integers M-N (required, scientific notation supported)\n"
        "       t             integer t to use for single N (optional, a good t will be determined if unspecified)\n");

}

int main (int argc, char *argv[])
{

    if ( argc < 2 ) { usage(); return 0; }

    int verbosity=0, optimal=0, create=0;
    char *hintfile = 0;
    int64_t minN=0, maxN=0, t=0;

    for ( int i = 1 ; i < argc ; i++ ) {
        char *s = argv[i];
        if ( maxN > minN || t ) { fprintf(stderr, "ignoring extraneous argument %s\n",s); continue; }
        if ( *s == '-' ) {
            switch(*(s+1)) {
            case 'v': { if ( i+1 >= argc) { usage(); return -1; } verbosity = atoi(argv[i+1]); i++; break; }
            case 'h': { if ( i+1 >= argc || hintfile ) { usage(); return -1; } hintfile = argv[i+1]; i++; break; }
            case 'c': { create = 1; break; }
            case 'o': { optimal = 1; break; }
            default: { usage(); return -1; }
            }
        } else {
            if ( !minN ) {
                if (*s=='[')s++;
                double x = strtod(s,&s);
                if ( (x-round(x)) > 0.0001 ) { fprintf(stderr, "N=%f must be an integer.\n",x); usage(); return -1; }
                minN = round(x);
                if ( s && *s ) {
                    if ( *s == '.' ) { while (*(++s)=='.'); }
                    else if ( *s == '-' || *s == ',' ) { s++; }
                    else { puts(s); usage(); return -1; }
                    maxN = strtold(s,0);
                    assert(maxN >= minN);
                } else {
                    maxN = minN;
                }
            } else {
                t = strtold(s,0);
            }
        }
    }

    if ( minN < 10 || maxN >= MAXN ) { fprintf(stderr,"N-range [%ld,%ld] must be contained in [10,2^40)\n", minN, maxN); return -1; }
    if ( t && t < cdiv(minN,3) ) { fprintf(stderr,"t=%ld must be at least N/3\n", t); return -1; }


    double start = get_time();
    general_setup();
    if ( verbosity > 0 ) fprintf(stderr,"General setup took %.3fs\n", get_time()-start);

    start = get_time();
    if ( maxN > minN ) {
        if ( ! hintfile ) { fprintf(stderr, "You must specify a hint-file (either to read or to create) when using a range of N.\n"); return -1; }
        if ( create ) {
            if ( !hintfile ) { fprintf(stderr, "You must use the -h parameter to specify the hint-file to be created.\n"); return -1; }
            FILE *fp = fopen(hintfile, "w");
            if ( !fp ) { fprintf(stderr, "Error creating hint-file %s\n", hintfile); return -1; }
            int64_t N = minN;
            while ( N <= maxN ) {
                t = tbound(N,optimal,verbosity);
                if ( !t ) break;
                if ( verbosity >= 0 ) printf ("t(%ld) >= %ld (t-N/3 >= %ld) (%.3fs)\n", N, t, t-cdiv(N,3), get_time()-start);
                fprintf(fp,"%ld:%ld\n",N,t);
                N += t-cdiv(N,3)+1;
            }
            if ( N > maxN ) fprintf (stdout,"Verified the Guy-Selfridge-Erdos conjecture for N in [%ld,%ld] (%.3fs)\n",minN,maxN,get_time()-start);
            else if ( N == minN ) fprintf (stdout,"Unable to verify Guy-Selfridge-Erdos conjecture for N=%ld (%.3fs)\n",minN,get_time()-start);
            else fprintf (stdout,"Only able to verify the Guy-Selfridge-Erdos conjecture for N in [%ld,%ld] (%.3fs)\n",minN,N-1,get_time()-start);
        } else {
            FILE *fp = fopen(hintfile,"r"); if (!fp) { fprintf(stderr, "Error opening hint-file %s\n", hintfile); return -1; }
            char buf[256];
            int64_t minV = 0, maxV = 0;
            while ( fgets(buf,sizeof(buf),fp) ) {
                int64_t N = atol(buf);
                char *s = strchr(buf,':'); if (!s) { fprintf(stderr, "Error parsing line %s\n", buf); return -1; }
                t = atol(s+1);
                if ( 3*t < N ) { fprintf(stderr, "Invalid N:t in hint file: 3*%ld < %ld\n", t, N); return -1; }
                double timer = get_time();
                if ( tfac(N,t,0,verbosity) < N ) { fprintf (stderr, "Failed to verify t(%ld) >= %ld !\n", N,t); return -1; }
                if (!minV) {
                    if ( N > minN ) { fprintf(stderr, "Hint file starting N=%ld above range minimum %ld\n", N, minN); return -1; }
                    minV = N;
                    maxV = N + (t-cdiv(N,3));
                } else {
                    if ( N > maxV+1 ) { fprintf(stderr, "Hint file starting N=%ld leaves a gap!\n", N); return -1; }
                    if ( N + (t-cdiv(N,3)) <= maxV ) { fprintf (stderr, "Hint at N=%ld did not extend verified range!\n", N); return -1; }
                    maxV = N + (t-cdiv(N,3));
                }
                if ( verbosity >= 0 ) printf ("t(%ld) >= %ld (%.3fs)\n", N, t, get_time()-timer);
                if ( maxV >= maxN ) break;
            }
            if ( maxV < maxN ) { fprintf (stderr, "Hint file only allowed verification [%ld,%ld]\n", minV,maxV); return -1; }
            fprintf (stdout,"Verified the Guy-Selfridge-Erdos conjecture for N in [%ld,%ld] (%.3fs)\n",minV,maxV,get_time()-start);
        }
    } else {
        int64_t N = minN;
        if ( t && optimal ) { t=0; fprintf(stderr,"Ignoring specified value of t and searching for optimal value\n"); }
        if ( !t ) {
            t = tbound(N,optimal,verbosity);
            if ( t ) printf("t(%ld) >= %ld%s with (t-ceil(N/3)) = %ld (%.3fs)\n",N,t,optimal ? " (optimal for algorithm)" : "", (t-cdiv(N,3)), get_time()-start);
            else fprintf(stderr,"failed to prove t(%ld) >= %ld (%.3fs)\n", N, t, get_time()-start);
        } else {
            int64_t cnt = tfac(N,t,0,verbosity);
            if ( cnt >= N ) printf("t(%ld) >= %ld with %ld extra factors (%.3fs)\n", N, t, cnt - N, get_time()-start);
            else fprintf (stderr,"failed to prove t(%ld) >= %ld with %ld missing factors (%.3fs)\n", N, t, N - cnt, get_time()-start);
        }
    }

    return 0;
}
