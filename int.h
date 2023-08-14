#ifndef _INT
#define _INT

// Integer utility routines and type definitions

#include<cstdint>
#include<cmath>
#include<cstdint>

#include "bigint.h"


using namespace std;

// types

typedef int16_t int16;
typedef int32_t int32;
typedef uint64_t uint64;
typedef int64_t int64;
typedef char bit;
// we have int128 from bigint.h


const int primeslen = 168;
const int16 primesmax = 1000;
const int16 primes[] = {
   2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
  73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
 661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
 947, 953, 967, 971, 977, 983, 991, 997
};

// testing for squares
bool issquare(int64 x)
{
  int64 root = (int64) llroundl(sqrtl((long double)x));
  return (root*root == x);
}

// GCD and inverse functions
int64 gcd(int64 x, int64 y)
{
  if(x<0) x= -x; if(y<0) y= -y;
  int64 r;
  while(y!=0) { r=x%y; x=y; y=r; }
  return x;
}

inline int64 lcm(int64 a, int64 b)
  { return (a/gcd(a,b))*b; }

int64 extgcd(int64 a, int64 b, int64 &x, int64 &y)
{
  int64 ux=1,uy=0,vx=0,vy=1,u=a,v=b,r,rx,ry,q;
  int64 asign=1, bsign=1;
  if(a<0) { a=-a; asign=-1; }
  if(b<0) { b=-b; bsign=-1; }
  while(v>0)
  {
    q=u/v; r=u-q*v;
    rx=ux-q*vx; ux=vx; vx=rx;  ry=uy-q*vy; uy=vy; vy=ry; 
    u=v; v=r;
  }
  x=asign*ux; y=bsign*uy;
  return u;
}

// returns the inverse of x modulo m
int64 inv(int64 x, int64 m)
{
  int64 a,b,g;
  g=extgcd(x%m,m,a,b);
  if(g!=1) return 0;
  if(a>0 && a<m) return a;
  if(a>0) return a%m;
  if(a>-m) return m+a;
  return m-((-a)%m);
}

// modular exponentiation: 32 and 64-bit versions

int32 powmod(int32 xin, int32 e, int32 m)
{
  int64 y=1;
  int64 x=xin;
  while(e>0)
  {
    if(e%2==1) y=(y*x)%m; // if e is odd
    e >>= 1;  // e = e/2
    x = (x*x)%m;
  }
  return (int32)y;
}

int64 powmod(int64 xin, int64 e, int64 m)
{
  int128 y=1;
  int128 x=xin;
  while(e>0)
  {
    if(e%2==1) y=(y*x)%m; // if e is odd
    e >>= 1;  // e = e/2
    x = (x*x)%m;
  }
  return (int64)y;
}

inline int16 legendre(int32 a, int32 p)
  { 
	int16 l=powmod(a, (p-1)/2, p);
	if (l==p-1)
		l=-1;
	return l; 
  }
inline int16 legendre(int64 a, int64 p)
  { 
	int16 l=powmod(a, (p-1)/2, p);
	if (l==p-1)
		l=-1;
	return l; 
  }

 inline int16 jacobi(int64 a, int64 b)
{
    int64 t=1;
    while (a!=0)
    {
        while ((a%2)==0)
        {
            a=a/2;
            if((b%8)==3 || (b%8)==5)
            {
                t=-t;
            }
        }
        swap(a, b);
        if((a%4)==3 && (b%4)==3)
        {
            t=-t;
        }
        a=(a%b); 
    } 
    if(b==1) {return t;}
    else {return 0;}
}


void createExpArr(int64 *Exponents, int64 *Primes, int64 Primelen, int64 num)
{
  int n=num;

  for(int len=0; len<Primelen; len++)
  {
    while((n%Primes[len])==0)
    {
      Exponents[len]+=1;
      n = n / Primes[len];
    }
  }
}

void Merge(int64 *p_dP, int64 *p_PD, int64 *e_dP, int64 *e_PD, int64 *prim_merge, int64 *exp_merge, int64 primes_dP_len, int64 primes_PD_len, int64& len)
  {
    int i=0;
    int j=0;
    int k=0;
    len=0;


    while(i<primes_dP_len && j<primes_PD_len)
    {
        if(p_dP[i] < p_PD[j])
        {
            if(e_dP[i]>0)
            {
                prim_merge[len] = p_dP[i];
                exp_merge[len] = e_dP[i];
                len++;
            }
            i++;
        }
        else if(p_dP[i] == p_PD[j])
        {
            int new_exp = e_dP[i] + e_PD[j];
            if( new_exp > 0 )
            {
                prim_merge[len] = p_dP[i];
                exp_merge[len] = new_exp;
                len++;
            }
            i++;
            j++;
        }
        else
        {
            if(e_PD[j]>0)
            {
                prim_merge[len] = p_PD[j];
                exp_merge[len] = e_PD[j];
                len++;

            }
            j++;
        }
    }

    while(i<primes_dP_len)
    {
        if(e_dP[i]>0)
        {
            prim_merge[len] = p_dP[i];
            exp_merge[len] = e_dP[i];
            len++;
        }
        i++;
    }

    while(j<primes_PD_len)
    {
        if(e_PD[j]>0)
        {
            prim_merge[len] = p_PD[j];
            exp_merge[len] = e_PD[j];
            len++;
        }
        j++;
    }

}


void getDivisors(int64 *primes, int64 *exponents, int64 initial_div, int64 primeslen, int64 *divisors, int64& z)
  {
    int64 div_count = 1;
    for(int y = 0; y < primeslen; y++)
    { 
        div_count*=(exponents[y]+1);
    }

    divisors[div_count];
    divisors[0]=initial_div;

    int64 currentloc; 
    int64 currdiv_count=1;
    z=1;
    for(int i=0; i<primeslen; i++)
    {
        int64 multiplier=1;
        for(int j=0; j<exponents[i]; j++)
        {
            multiplier*=primes[i];
            for(int k=0; k<currdiv_count; k++)
            {
                divisors[z]=  multiplier * divisors[k];
                z++;
            }
        }
        currdiv_count*=(exponents[i]+1);
    }
    
  } 

bool KorseltCrit(int64 *current, int64 currentval, int64 length, int64 q, int64 eq, int64 er, int64 r, int64 eP)
  {
    int64 en = eP*eq*er;

    int64 temp_residueQ;
    int64 temp_residueR;
    int64 temp_residueP;

    temp_residueQ =  (currentval * r) % ( q - eq );
    temp_residueQ = (temp_residueQ * q) % ( q - eq );
    temp_residueQ = (temp_residueQ - en) % ( q - eq );

  
    if( temp_residueQ !=0 )
    {
      return false;
    }
    
    temp_residueR = (currentval * r) % (r - er);
    temp_residueR = (temp_residueR * q) % (r - er);
    temp_residueR = (temp_residueR - en) % (r - er);

    if( temp_residueR !=0 )
    {
      return false; 
    }
    
    for(int i=0; i<length; i++)
    {
        int ec = -1;
        if(current[i]%5==1 || current[i]%5==4)
        {
            ec = 1;
        }

        temp_residueP = (currentval * r) % (current[i] - ec);
        temp_residueP = (temp_residueP * q) % (current[i] - ec);
        temp_residueP = (temp_residueP - en) % (current[i] - ec);
        if( temp_residueP !=0 )
        {
            return false;
        
        }
    }
    
    return true;
   
  }
  bool KorseltCrit2(int64 *Primes_dividing_P, int64 *Sign, int64 length, uint64 ndn)
  {
    for(int i=0; i<length; i++)
    {
      if(ndn % (Primes_dividing_P[i]-Sign[i]) != 0)
      {
        return false;
      }
      
    }
    return true;

  }



  bool sprp(int64 n, int64 a) {
  /* Calculate d/s representation of n */
  int64 d=n-1;
  int16 s=0;
  while (!(d & 0xff)) { d>>=8; s+=8; }
  if (!(d & 0xf)) { d>>=4; s+=4; }
  if (!(d & 0x3)) { d>>=2; s+=2; }
  if (!(d & 0x1)) { d>>=1; s+=1; }
  // Calculate a^d(mod n)
  int64 b=powmod(a,d,n);
  if ((b==1) || (b==(n-1))) return true;
  int16 r;
  for (r=1; r<s; r++) {
    b=b*b % n;
    if (b<=1) return false;
    if (b==(n-1)) return true;
  }
  return false;
}


bool isPrime(uint64 n) {
  // Catch easy answers
  if (n<2) return false;    // 0 and 1 are not prime
  if (n<4) return true;     // 2 and 3 are prime
  if (!(n&1)) return false; // Even numbers are not

 
  for (int i=0; i<primeslen; i++) {
    int64 p_test = primes[i];
    if (n==p_test) return true;
    if (n % p_test == 0) return false;
  }

  // Next step, SPRP tests
  //
  // Thresholds from Sloan sequence A014233
  if (!sprp(n,2)) return false;
  //if (n<2047) return true;
  if (!sprp(n,3)) return false;
  //if (n<1373653) return true;
  if (!sprp(n,5)) return false;
  if (n<25326001) return true;
  if (!sprp(n,7)) return false;
  //if (n<3215031751ULL) return true;
  if (n==3215031751) return false;
  if (n<118670087467ULL) return true;
  if (!sprp(n,11)) return false;
  if (n<2152302898747ULL) return true;
  if (!sprp(n,13)) return false;
  if (n<3474749660383ULL) return true;
  if (!sprp(n,17)) return false;
  if (n<341550071728321ULL) return true;
  if (!sprp(n,19)) return false;
  if (!sprp(n,23)) return false;
  if (n<3825123056546413051ULL) return true;
  if (!sprp(n,29)) return false;
  if (!sprp(n,31)) return false;
  if (!sprp(n,37)) return false;
  // This test passes for n<2^64
  return true;
}


void Prime_init(int64 primes_vec, int64 prime_symbol, int64 *Primes_dividing_P, int64 *Prime_symbol, int64 *Sign, int64 *Running_L, int64 *Running_Prod, int64 *Running_Sign)
{
      Primes_dividing_P[0]=primes_vec;
      Prime_symbol[0]=primes_vec-prime_symbol;
      Sign[0]=prime_symbol;
      Running_L[0]=Prime_symbol[0];
      Running_Prod[0]=primes_vec;
      Running_Sign[0]=Sign[0];
}

inline void Prime_update(int64 primes_vec, int prime_symbol, int64 v, int64 *Primes_dividing_P, int64 *Prime_symbol, int64 *Sign, int64 *Running_L, int64 *Running_Prod, int64 *Running_Sign, int64 prev)
    {
      Primes_dividing_P[v]=primes_vec;
      Prime_symbol[v]=primes_vec-prime_symbol;
      Sign[v]=prime_symbol;
      Running_L[v]=lcm(Prime_symbol[v], Running_L[prev]);
      Running_Prod[v]=primes_vec*Running_Prod[prev];
      Running_Sign[v]=Sign[v]*Running_Sign[v-1];
    }


#endif
