/*
  Write N! as the product of N integer factors.  How large can the
  smallest factor be?

  This program generates a linear program to help solve the problem.
  The output is in the "LP file format", or "lp-format", which is
  lpsolve's native format.

  Note that the linear program generated assumes the factors are all
  at most N, which has not (yet) rigorously shown to be a valid
  assumption.  (In fact, for T larger than N/2, it is not valid.)

  Compile with something like:
  g++ -O2 -std=c++14 generate_lp_format.cc -o generate_lp_format

  Run with something like:
  ./generate_lp_format -i 7 2 | lp_solve

  https://arxiv.org/abs/2503.20170
  https://terrytao.wordpress.com/2025/03/26/decomposing-a-factorial-into-large-factors/
  https://en.wikipedia.org/wiki/Lp_solve
 */

#include <iostream>
#include <cstring>
#include <vector>

int main( int argc, char **argv ) {
	// Print usage message
	argc--; argv++;
	if( argc == 0 or std::strchr(*argv, 'h') != nullptr ) {
		std::cerr << "Usage: generate_lp_format [-f] [-i] N [T]\n\n";
		std::cerr << "N! is the factorial to be factored into >=N factors.\n";
		std::cerr << ">=T is the threshold for how large the factors must be, ceil(N/3) by default.\n";
		std::cerr << "-f merely tests feasibility rather than maximizing the number of factors.\n";
		std::cerr << "-i enforces that the decision variables are integers.\n";
		return 1;
	}

	// Parse arguments, very/too simply
	bool feasibility_only = false;
	if( argc > 0 and std::strchr(*argv, 'f') != nullptr ) {
		argc--; argv++;
		feasibility_only = true;
	}

	bool integral_variables = false;
	if( argc > 0 and std::strchr(*argv, 'i') != nullptr ) {
		argc--; argv++;
		integral_variables = true;
	}

	int N = 0;
	if( argc > 0 ) {
		N = atoi(*argv);
		argc--; argv++;
	}
	int T = (N+2)/3;
	if( argc > 0 ) {
		T = atoi(*argv);
		argc--; argv++;
	}
	if( not(N > 0 and T > 0 and N >= 2*T) ) {
		std::cerr << "Unexpected N and T.  Maybe command-line argument interpreted incorrectly.\n";
		exit(10);
	}

	// Compute primes up to N with a simple sieve of Eratosthenes
	std::vector<bool> is_prime(N+1, true);
	is_prime[0] = is_prime[1] = false;
	std::vector<int> primes;
	for( int p = 2; p <= N; p++ ) {
		if( not is_prime[p] )
			continue;
		primes.push_back(p);
		for( int i = 2*p; i <= N; i += p ) {
			is_prime[i] = false;
		}
	}

	if( feasibility_only ) {
		// no objective
		std::cout << "max: ;\n";
		// ... but we do insist there be N factors
		for( int x = T; x <= N; x++ ) {
			std::cout << "+x" << x;
		}
		std::cout << " >= " << N << ";\n";
	} else {
		std::cout << "max: ";
		for( int x = T; x <= N; x++ ) {
			std::cout << "+x" << x;
		}
		std::cout << ";\n";
	}

	for( int p : primes ) {
		std::cout << "c" << p << ": ";
		for( int x = T; x <= N; x++ ) {
			// Compute the valuation of x at p
			int v = 0;
			int t = x;
			while( (t % p) == 0 ) {
				v++;
				t /= p;
			}
			if( v > 0 ) {
				std::cout << "+";
				if( v > 1 )
					std::cout << v;
				std::cout << "x" << x;
			}
		}
		int rhs = 0;
		{
			int t = N;
			while( t ) {
				t /= p;
				rhs += t;
			}
		}
		std::cout << " <= ";
		std::cout << rhs;
		std::cout << ";\n";
	}

	if( integral_variables ) {
		std::cout << "int ";
		for( int x = T; x <= N; x++ ) {
			if( x > T )
				std::cout << ",";
			std::cout << "x" << x;
		}
		std::cout << ";\n";
	}
}
