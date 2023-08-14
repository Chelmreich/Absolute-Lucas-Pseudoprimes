#ifndef _PT
#define _PT

// functions for deterministic, unconditional prime testing
// for small-ish integers

#include "int.h"

using namespace std;


// 32 and 64-bit strong pseudoprime test for the given base
// we have n-1=m*2^e with m odd
bool strong(int32 base, int32 n, int32 m, int16 e)
{
  int64 x=powmod(base, m, n);
  if(x==1 || x==n-1) return true;
  for(int i=0; i<e; i++)
  {
    x=(x*x)%n;
    if(x==n-1) return true;
    if(x==1) return false;
  }
  return false;
}
bool strong(int64 base, int64 n, int64 m, int16 e)
{
  int128 x=powmod(base, m, n);
  if(x==1 || x==n-1) return true;
  for(int i=0; i<e; i++)
  {
    x=(x*x)%n;
    if(x==n-1) return true;
    if(x==1) return false;
  }
  return false;
}

// fast test for 16-bit integers
// ψ2 = 1373653 > 2^15=32768
// so two strong ps tests will do it
bool primetest(int16 n)
{
//cout << "Prime test 16 bit of " << n << endl;
  if(n==1) return false;
  if(n==2 || n==3 || n==5) return true;
  if(n%2==0) return false; // n is even
  int16 m=(n-1)/2;
  int16 e=1;
  while(m%2==0) { e++; m=m/2; }
  if(!strong(2,n,m,e)) return false;
  return strong(3,n,m,e);
}

// fast test for 32-bit integers
// ψ4 = 3215031751 > 2^31=2147483648
// so four strong ps tests will do it
bool primetest(int32 n)
{
//cout << "Prime test 32 bit of " << n << endl;
  if(n==1) return false;
  if(n==2 || n==3 || n==5 || n==7) return true;
  if(n%2==0) return false; // n is even
  int32 m=(n-1)/2;
  int16 e=1;
  while(m%2==0) { e++; m=m/2; }
  if(!strong(2,n,m,e)) return false;
  if(!strong(3,n,m,e)) return false;
  if(!strong(5,n,m,e)) return false;
  return strong(7,n,m,e);
}

const int64 primetestcutoff=3825123056546413051;

// test for 64-bit integers
// ψ9 = 3825123056546413051, below 2^63=9223372036854775808,
// so we punt above this limit, for now
bool primetest(int64 n)
{
//cout << "Prime test 64 bit of " << n << endl;
  if(n<2147483648) return primetest((int32) n);
  if(n%2==0) return false; // n is even
  if(n<primetestcutoff) // we use 9 strong ps tests
  {
    int64 m=(n-1)/2;
    int16 e=1;
    while(m%2==0) { e++; m=m/2; }
    if(!strong((int64) 2,n,m,e)) return false;
    if(!strong((int64) 3,n,m,e)) return false;
    if(!strong((int64) 5,n,m,e)) return false;
    if(!strong((int64) 7,n,m,e)) return false;
    if(!strong((int64) 11,n,m,e)) return false;
    if(!strong((int64) 13,n,m,e)) return false;
    if(!strong((int64) 17,n,m,e)) return false;
    if(!strong((int64) 19,n,m,e)) return false;
    return strong((int64) 23,n,m,e);
  }
  else 
    return false; // for now...
}

#endif
