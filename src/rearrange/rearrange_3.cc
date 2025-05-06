/*
This program implements a dynamic programming approach to rearranging
powers of 2 and 3 in the standard factorization of N!.  On
command-line input "PROGRAM N T", it prints whether or not it is
possible to get t(N) >= T by rearranging only the tiny primes 2,3.
(In other notation, this tests whether t_{2,3}(N) >= T.)

More or less, after removing all powers of 2 and 3, it computes the
optimal tradeoff curve between powers of 2 and 3 in terms of bringing
up the rest of the factors up to T.
 */

#include <algorithm>
#include <iostream>
#include <vector>

bool possible( int N, int T ) {
	// Track the amount of powers of 2 and 3 available
	int budget2 = 0;
	int budget3 = 0;
	std::vector<int> target(0);
	for( int i = 1; i <= N; i++ ) {
		int t = i;
		while( (t % 2) == 0 ) {
			t /= 2;
			budget2++;
		}
		while( (t % 3) == 0 ) {
			t /= 3;
			budget3++;
		}
		// (T+t-1)/t is ceil(T/t), the smallest integer we can
		// multiply t by to get T.
		target.push_back((T+t-1) / t);
	}
	std::sort( target.begin(), target.end() );

	// Dynamic program to determine whether our budgets suffice:
	// two[three] counts the smallest power of 2 that suffices for
	// a given budget of powers of 3, in terms of multiplying up
	// all numbers considered so far.
	std::vector<int> two(budget3+1, 0);

	// here_two[here_three] does the same, but only for the
	// current target number t.
	int max3 = 0;
	for( int pow3 = 1; pow3 <= N; pow3 *= 3 ) {
		max3++;
	}
	std::vector<int> here_two(max3+1, 0);
	std::vector<int> helper  (max3+1, 1);
	for( int j = 1; j <= max3; j++ ) {
		helper[j] = helper[j-1]*3;
	}
	for( int t : target ) {
		if( t == 1 )
			continue;
		// We are trying to get a factor of at least t.  How
		// many powers do we need here?
		for( int j = 0; j <= max3; j++ ) {
			while( helper[j] < t ) {
				here_two[j]++;
				helper[j] *= 2;
			}
		}
		// Main recursive step of the dynamic program
		for( int j = budget3; j >= 0; j-- ) {
			two[j] += here_two[0];
			for( int k = 1; k <= max3 and j-k >= 0; k++ ) {
				if( two[j] > two[j-k] + here_two[k] ) {
					two[j] = two[j-k] + here_two[k];
				}
			}
		}
	}

	return (two[budget3] <= budget2);
}

int main(int argc, char **argv) {
	int N = 0;
	int T = 0;
	argc--; argv++;
	if( argc ) {
		N = atoi(*argv);
		argc--; argv++;
	}
	if( argc ) {
		T = atoi(*argv);
		argc--; argv++;
	}
	if( not(N > 0 and T > 0) ) {
		std::cerr << "Usage: PROGRAM N T\n";
		return 1;
	}

	bool result = possible(N, T);
	std::cout << int(result) << std::endl;

	return 0;
}
