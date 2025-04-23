#include <iostream>
#include <vector>
using namespace std;

// these constants define the interval for which we run the program
// we require minN*minN >= maxN
#define minN 10000000LL
#define maxN 100000000000LL

// this is the number of consecutive intervals that are used for precomputing primes and the factorization
#define NUMBEROFINTERVALS 10000

// this is the largest itneger for which we compute factorization (for larger integer we factorize it only if it is divisible by 2 or 3)
#define LARGEST_COMPUTED_FACTORIZATION maxN/200

vector<long long> primes;
int* smallestDivisor = new int[LARGEST_COMPUTED_FACTORIZATION/2+1];

///////////////////////////////////
//// Precomputations of primes ////
///////////////////////////////////

//compute the primes numbers and factorization for interval [0, e]
void factorizationOfFirstInterval(long long e){
  vector<long long> firstDivisor(e+1, -1);
  for(long long i=2; i<firstDivisor.size(); i++){
    if(firstDivisor[i] == -1){
      primes.push_back(i);
      firstDivisor[i] =  primes.size()-1;
      if(i <= LARGEST_COMPUTED_FACTORIZATION && i%2 == 1){
        smallestDivisor[i/2] = firstDivisor[i];
      }
      for(long long j = i; j*i<firstDivisor.size(); j++){
        if(firstDivisor[j*i] == -1){
          firstDivisor[j*i] = primes.size()-1;
          if(j*i <= LARGEST_COMPUTED_FACTORIZATION && (i*j)%2 == 1){
            smallestDivisor[(j*i)/2] = primes.size()-1;
          }
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
      if(firstDivisor[j*p-b] == -1){
        firstDivisor[j*p-b] = i;
      } 
    } 
  }
  for(long long i=b; i<=e; i++){
    if(i%1000000000 == 0) cout <<"Precomputation reached  i ="<< i << " out of " << maxN << endl;
    if(firstDivisor[i-b] != -1){
        long long fd = primes[firstDivisor[i-b]];
        if( i > LARGEST_COMPUTED_FACTORIZATION || i%2 == 0) continue;         
        smallestDivisor[i/2] = firstDivisor[i-b];
    } else{
      primes.push_back(i);
      if( i > LARGEST_COMPUTED_FACTORIZATION || i%2 == 0) continue;           
      smallestDivisor[i/2] = primes.size()-1;
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
      for (size_t i = 0; i < /*factors.size()*/ 100000; ++i) {
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
      //this is quit ugly, but this avoids overflow in the multiplication
      //other solution would involve anoying while conditions
      __int128 d = p;
      while(d <= n){
        factors[i] += n/d;
        d=d*p;
      }
    }
  }
   
  //tbr is the index of one prime to be removed
  //toberemoved is a number to be removed (so we need to take its factorization)
  long long countRemovable(long long toBeRemoved, long long tbr){
    if(toBeRemoved%2 != 0 && toBeRemoved%3 != 0 && toBeRemoved > LARGEST_COMPUTED_FACTORIZATION ) return 0;
    // we remove the content of the vector and the integer tbr
    // we use the fact that in toBeRemoved all occurences of the same element are next to each other
    long long removable = factors[tbr];

    while(toBeRemoved>1){
      long long v;
      if(toBeRemoved%2 == 0) v = 0;
      else if(toBeRemoved%3 == 0) v = 1;
      else if(toBeRemoved > LARGEST_COMPUTED_FACTORIZATION) return 0;
      else v = smallestDivisor[toBeRemoved/2];
      toBeRemoved /= primes[v];      
      long long nbocc=1;
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
    found += nbTimes;
    factors[tbr] -= nbTimes;
    long long val = toBeRemoved*primes[tbr];
    if(val < targetVal){    //sanity check 1
      cout<<"error too small"<<endl;
      exit(1);
    }
    while(toBeRemoved>1){
      long long v;
      if(toBeRemoved%2 == 0) v = 0;
      else if(toBeRemoved%3 == 0) v = 1;
      else v = smallestDivisor[toBeRemoved/2];
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
    while(countAvailable == 0 && goal < maxN/20){
      goal ++;
      countAvailable = factorizer.countRemovable(goal, i);
    }
    prevGoal = goal;
    if(goal >= maxN/20) break;
    factorizer.addToFactorization(countAvailable, goal, i);
  }
  return factorizer.found;
}


int main(){
  if(minN < maxN/minN){
    cout<<"Error: minN*minN < maxN"<<endl;
    return 0;
  }
  cout<<"Running with range [ "<<minN<<" , "<<maxN<<" ]"<<endl;
  cout<<"Precomputing the primes and factorizations."<<endl;
  setUpFactorization();
  cout<<"There are "<<primes.size() <<" primes lessor equal than "<<maxN<<endl;
  factorizer.setMemory();
  vector<long long> eps = {30,20,15,10,5,4,3,2,1,0};
  for(long long n=minN; n<=maxN;){
    if(n%3 != 0 && n>minN+1){ n++; continue;}
    bool done = false;
    for(long long e=0; e<eps.size() && !done; e++){
      if(n<10000000000LL && eps[e]>40) continue;
      long long res = bestFact(n, eps[e]);
      if(res >= n){
        done = true;
        cout << "Donne for n = "<<n<<" with eps = "<<eps[e]<< " next n = "<<(n*(1000+eps[e])/1000)+1<<endl;
        cout << "Each of the factors found where larger than: " << factorizer.targetVal << " (n/3 = "<< (n/3) <<")."<< endl;
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
