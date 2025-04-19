#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

//the following might not be available on your machine (it should on a linux machine whit compiler g++)
//in this case simply replace the cc_hash_table by unordered_map which are slower
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
cc_hash_table<long long, long long> smallestDivisor;
cc_hash_table<long long, long long> largestDivisor;

// unordered_map<long long, long long> smallestDivisor;
// unordered_map<long long, long long> largestDivisor;


// these constant define the interval for which we run the program
// we require minN*minN >= maxN
#define minN 10000000LL
#define maxN 100000000000LL

// this is the number of consecutive intervals that are used for precomputing primes and the factorization
#define NUMBEROFINTERVALS 400


// Chosen so that the product of the entries in a given list is not too large.
// The current choice seems to work well. 
//
// Each list corresponds to a single array lookup for the primes covered.
// So, e.g. checking divisibilty one prime at a time would roughly correspond to having
// a separate list for each prime. For convenience, here is the % coverage for including all primes
// up to p for the first several values of p:
//  2 : 50.0 %
//  3 : 66.7 %
//  5 : 73.3 %
//  7 : 77.1 %
//  11 : 79.2 %
//  13 : 80.8 %
//  17 : 81.9 %
//  19 : 82.9 %
//  23 : 83.6 %
//  29 : 84.2 %
// 
//  Note that getting above 90% requires the first 55 primes.
// 
// In general, more lists means more run-time lookups, but less memory used and hence better cacheability.
// These can be rearranged or edited for tradeoffs. As written, the only logical requirement is they must 
// cover exactly the first N primes for some N, and the primes must appear in order.
const vector<vector<long long>> small_primes{
    {2, 3, 5, 7}, {11, 13} 
    // {17, 19},
    // {23, 29, 31, 37},
    // {41, 43, 47, 53}
}; 

// After running prepSmallPrimes, this will contain the same number of vectors as 
// small_primes, but each one has length equal to the product of the primes in the corresponding
// list in small_primes. The ith term of a sub-vector of small_prime_mods is
// equal to the index in primes of the smallest prime that both divides i and belongs to the 
// corresponding list in small_primes. Otherwise, -1 if i is not divisible by any of those primes. 
// Thus, for any integer i > 0, if the smallest divisor of i is one of the small_primes,
// the value can be determined by looking up, in order, the value of v[i%v.size()] for each v in small_prime_mods.
vector<vector<long long>> small_prime_mods;

vector<long long> primes;

void prepSmallPrimes(){
  // Populate directly as we will skip small primes in the main loop. 
  for (const auto& v : small_primes) {
    for (auto p : v) {
      primes.push_back(p);
    }
  }
  
  long long curr_index = 0;

  for (int i = 0; i < small_primes.size(); ++i){
    long long vsize = 1;
    const auto& curr_prime_list = small_primes[i];
    // Compute the product of the primes in the sublist.
    for (int j = 0; j < curr_prime_list.size(); ++j) {
      vsize *= curr_prime_list[j];
    }
    // Allocate the vector of that size
    vector<long long> curr_prime_mod(vsize, -1);
    // Populate the vector
    for (int j = 0; j < curr_prime_list.size(); ++j) {
      auto p = curr_prime_list[j];
      for(int k = 0; p*k<vsize; ++k) {
        if(curr_prime_mod[p*k] == -1) {
          curr_prime_mod[p*k] = curr_index;
        }
      }
      curr_index++;
    }

    small_prime_mods.push_back(curr_prime_mod);

  }
  if(curr_index != primes.size()) {
    cout << "error in small prime indexing" << endl;
    exit(1);
  }
}

// Returns the index into primes of the smallest divisor of x.
long long getSmallestDivisor(long long x){
  for (const auto& v : small_prime_mods) {
    auto val = v[x%v.size()];
    if (val != -1) {
      return val;
    }
  }
  return smallestDivisor[x]; 
}

bool hasSmallDivisor(long long x){
  for (const auto& v : small_prime_mods) {
    auto val = v[x%v.size()];
    if (val != -1) {
      return true;
    }
  }
  return false; 
}

///////////////////////////////////
//// Precomputations of primes ////
///////////////////////////////////

//compute the primes numbers and factorization for interval [0, e]
//
// Note: We only populate smallestDivisor[i] and largestDivisor[i] 
// if i is not divisible by a small prime.
void factorizationOfFirstInterval(long long e){
  vector<long long> firstDivisor(e+1, -1);
  for(long long i=2; i<firstDivisor.size(); i++){
    if (hasSmallDivisor(i)){
      continue;
    }
    if(firstDivisor[i] != -1){
      long long fd = primes[firstDivisor[i]];
      // Note that i/fd cannot have a small divisor
      largestDivisor[i] = max(fd, largestDivisor[i/fd]);
    } else {
      primes.push_back(i);
      largestDivisor[i] = i;
      smallestDivisor[i] = firstDivisor[i] =  primes.size()-1;
      for(long long j = i; j*i<firstDivisor.size(); j++){
        if (hasSmallDivisor(j*i)) {
          continue;
        }
        if(firstDivisor[j*i] == -1){
          firstDivisor[j*i] = primes.size()-1;
          smallestDivisor[j*i] = primes.size()-1;
        } 
      }
    }
  } 
}

//compute the primes numbers and factorization for interval [b,e]
void factorizationOfOtherIntervals(long long b, long long e){
  vector<long long> firstDivisor(e+1-b, -1);
  for(long long i=0; primes[i]*primes[i]<=e; i++){
    long long p = primes[i];
    for(long long j = max(p, (b+p-1)/p); j*p<=e; j++){
      if (hasSmallDivisor(j*p)) continue;
      if(firstDivisor[j*p-b] == -1){
        firstDivisor[j*p-b] = i;
      } 
    } 
  }
  for(long long i=b; i<=e; i++){
    if(i%1000000000 == 0) cout <<"Precomputation eached  i ="<< i << "out of" << maxN << endl;
    if(hasSmallDivisor(i)) continue;
    if(firstDivisor[i-b] != -1){
        long long fd = primes[firstDivisor[i-b]];
        if(largestDivisor.find(i/fd) == largestDivisor.end()) continue;
        long long L = max(fd, largestDivisor[i/fd]);
        if(((i/fd)*3) >= maxN/L) continue;         
        smallestDivisor[i] = firstDivisor[i-b];
        largestDivisor[i] = L;
    } else{
      primes.push_back(i);
      if(i*3>maxN) continue;
      largestDivisor[i] = i;
      smallestDivisor[i] = primes.size()-1;
    }
  } 
}

// compute prime and factorization for the whole range by invoking the two previous functions
void setUpFactorization(){
  factorizationOfFirstInterval(minN);
  long long step = maxN/NUMBEROFINTERVALS;
  for(long long i=minN; i<maxN; i+=step){
    factorizationOfOtherIntervals(i+1, min(i+step, maxN));
  }
  largestDivisor.clear();
}



// a binary search for the index of the smallest prime larger or equal to v
// it is used once at the begining of every factorization of n!
// it is not clear to me whether it saves a lot of time
long long findIndexLargerPrime(long long v){
  long long u=primes.size()-1, l = 0;
  while(u-l>1){
    long long m = (u+l)/2;
    if(primes[m] > v) u = m;
    else l = m;
  }
  return u;
}


///////////////////////////////////
/////// Factorization of n! ///////
///////////////////////////////////

//A class that search the desired factorization of n!
//we really use it as a namespace more than a class
//currently we isntenciate only one object
//in particular the array factor is only created once (by the function setMemory)
class Factorizer{
  public:
  vector<long long> factors;
  long long targetVal;
  //the available factors from n!
  //the number of factors >targetVal already constructed
  long long  found=0;
  //list<pair<int, vector<long long>>> targetFactorization;
  //priority_queue<pair<long long,long long >, std::vector<pair<long long,long long >>, std::greater<pair<long long,long long >>> targetFactorization; 
  Factorizer(){
  };

  void printFactors(){
    bool atLeastOne = false;
    cout << "[";
      for (int i = 0; i < factors.size(); ++i) {
        if (factors[i] != 0) {
          if (atLeastOne) {
            cout << ", ";
          }
          atLeastOne = true;
          cout << primes[i] << ":" << factors[i];
        }
      }
    cout << "]" << endl;
  }

  //reserves enought memory for the array factors once and for all
  void setMemory(){
    factors.resize(primes.size());
  }
  //resets the array and the other attributes to get ready to factorize a new number
  void set(long long n, long long tv){
    targetVal = tv;
    found = 0;
    for(long long i=0; i<primes.size() && primes[i]<=n; i++){
      factors[i]=0;
      long long p = primes[i];
      long long n_tmp = n/p;
      while(n_tmp > 0){
        factors[i] += n_tmp;
        n_tmp/=p;
      }
    }
  }
  //tbr is the index of one prime to be removed
  //toberemoved is a number to be removed (so we need to take its factorization)
  long long countRemovable(long long toBeRemoved, long long tbr){
    if(!hasSmallDivisor(toBeRemoved) && smallestDivisor.find(toBeRemoved) == smallestDivisor.end()) return 0;
    // we remove the content of the vector and the integer tbr
    // we use the fact that in toBeRemoved all occurences of the same element are next to each other
    long long removable = factors[tbr];

    while(toBeRemoved>1){
      long long v;
      if(!hasSmallDivisor(toBeRemoved) && smallestDivisor.find(toBeRemoved) == smallestDivisor.end()) return 0;
      else v = getSmallestDivisor(toBeRemoved);
      toBeRemoved /= primes[v];      
      long long nbocc=1;
      // Note: The number of loop iterations here is v_p(toBeRemoved)
      // This could be modified to 1+log_2(v_p(toBeRemoved)) iterations by 
      // repeatedly squaring primes[v], then working your way back down. 
      while(toBeRemoved>1 &&  toBeRemoved % primes[v] == 0 ){
        nbocc++;
        toBeRemoved /= primes[v];
      }
      if(v == tbr) nbocc++;
      removable = min(removable, factors[v]/nbocc);
    }
    return removable;
  }
  // add the nbTimes the factorization of the number (toBeRemoved*primes[tbr])
  // this is always called with nbTimes large enough 
  // (most of the time after checking with a first call to countRemovable)
  // we still have two sanity checks, but they should not be necessary
  void addToFactorization(long long nbTimes, long long toBeRemoved, long long tbr){
    // Handling tbr part. 
    found += nbTimes;
    factors[tbr] -= nbTimes;
    long long val = toBeRemoved*primes[tbr];
    if(val < targetVal){    //sanity check 1
      cout<<"error too small"<<endl;
      exit(1);
    }
    // Handling toBeRemoved part
    // There might be some speedups here when toBeRemoved is divisible by a large prime power.
    while(toBeRemoved>1){
      long long v = getSmallestDivisor(toBeRemoved);
      toBeRemoved /= primes[v];      
      factors[v]-=nbTimes;
      if(factors[v]<0){    //sanity check 2
        cout<<"error neg val"<<endl;
        exit(1);
      }
    }
  }
  void printRes() const{
    cout<<"Number of factors found larger than "<<targetVal <<" is "<<found<<endl;
  }
};

Factorizer factorizer;

//this is the main function that use the previous class to factorize n! by using the greedy algorithm
long long bestFact(long long n, long long eps = 0){
  long long targetVal = (n*(1000+eps)/1000)/3;
  factorizer.set(n, targetVal);
  long long i = findIndexLargerPrime(n);
  for(; primes[i]>n; i--){}
  for(; primes[i]>=targetVal; i--){
    factorizer.addToFactorization(factorizer.factors[i], 1, i);
  }
  long long prevGoal = 2;
  while(i >= 0){
    if(factorizer.found >= n) return factorizer.found; 
    if(factorizer.factors[i] == 0){
      i--;
      continue;
    }

    long long p = primes[i]; 
    long long goal = max((targetVal+p-1)/p, prevGoal);  
    long long countAvailable = factorizer.countRemovable(goal, i);
    while(countAvailable == 0 && goal < targetVal){
      goal ++;
      countAvailable = factorizer.countRemovable(goal, i);
    }
    prevGoal = goal;
    if(goal >= targetVal) break;
    factorizer.addToFactorization(countAvailable, goal, i);
  }
  return factorizer.found;
}


int main(){
  prepSmallPrimes();
  if(minN < maxN/minN){
    cout<<"Error: minN*minN < maxN"<<endl;
    return 0;
  }
  cout<<"Running with range [ "<<minN<<" , "<<maxN<<" ]"<<endl;
  cout<<"Precomputing the primes and factorizations."<<endl;
  setUpFactorization();
  cout<<"There are "<<primes.size() <<" primes lessor equal than "<<maxN<<endl;
  factorizer.setMemory();
  vector<long long> eps = {60,50,40,30,20,15,10,5,4,3,2,1,0};
  for(long long n=minN; n<=maxN;){
    if(n%3 != 0 && n>minN+1){ n++; continue;}
    bool done = false;
    for(long long e=0; e<eps.size() && !done; e++){
      long long res = bestFact(n, eps[e]);
      if(res >= n){
        done = true;
        cout<<"Donne for n = "<<n<<" with eps = "<<eps[e]<< " next n = "<<(n*(1000+eps[e])/1000)+1<<endl;
        // Used to sanity check there is no difference in behavior from the previous version of the algorithm.
        cout << "Remaining factors: ";
        factorizer.printFactors();
        n = (n*(1000+eps[e])/1000)+1;
      }
    }
    if(!done){
      cout<<" !!!!   n = "<<n<<" NOT OK "<<endl;
      n++;
    }
  }
  return 0;
}
