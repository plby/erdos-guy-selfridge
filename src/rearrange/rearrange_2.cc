#include <iostream>
#include <queue>

int main( ) {
	int budget = 0;
	std::priority_queue<
		int,
		std::vector<int>,
		std::greater<int>
		> factors{};
	for( int N = 1; true; N++ ) {
		{
			int t = N;
			while( (t % 2) == 0 ) {
				t /= 2;
				budget++;
			}
			factors.push(t);
		}

		while( budget > 0 ) {
			budget--;
			int f = factors.top();
			factors.pop();
			factors.push(2*f);
		}
		int t = factors.top();
		std::cout << N << " " << t << std::endl;
	}
	
	return 0;
}
