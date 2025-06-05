#include <iostream>

int w( int N ) {
	return __builtin_popcount(N);
}

bool check( int N, int T ) {
	return (N >= 2 * ((8*T + 2 - 2*w(T)) / 3));
}

int main( ) {
	int T = 1;
	for( int N = 1; true; N++ ) {
		while( check(N,T) ) {
			T++;
		}
		std::cout << 2 << " " << N << " " << T << "\n";
	}

	return 0;
}
