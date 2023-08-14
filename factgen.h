#ifndef _FACTGEN
#define _FACTGEN

#include <cmath>
#include "int.h"
#include "primetest.h"
#include "stack.h"

using namespace std;



class Factgen
{

  Stack *roll;
  int rollsize;
  int pos;

public:

  Factgen(): roll(NULL) {}

  int64 prev[20];
  int prevlen;
  int64 prevn;
  int64 start, stop;
  int64 n;

  void print()
  {
    cout <<"n="<< n<<endl;
    cout <<"Empty?"<<roll[pos].isempty()<<endl;
    cout <<"stack top="<<roll[pos].gettop()<<endl;
  }

  void init(int64 startin, int64 stopin)
  {
    if(roll!=NULL) delete [] roll;
    rollsize=1+sqrt(stopin);
    rollsize=2*rollsize;
    roll=new Stack[rollsize];

    stop=stopin;
    start=startin;
    n=start;
    prevlen=0;
    prevn=0;
    pos=0;

    for(int32 p=rollsize-1; p>=2; p--)
      if(primetest(p))
      {
        roll[ (p-(start%p))%p ].push(p);
      }
//cout << "rollsize=" << rollsize << endl;
  }

  void next()
  {
    // copy stuff & factor
    int64 r=n;
    prevn=n; prevlen=0;
    while(!roll[pos].isempty())
    {
      int32 p=roll[pos].pop();
//cout << "pop " << p << endl;
      roll[(pos+p)%rollsize].push(p);
      prev[prevlen++]=p;
      while(r%p==0) { r=r/p; } // cout << "div\n"; }
    }
    if(r>1)
      { prev[prevlen++]=r; }
    n++;
    pos=(pos+1)%rollsize;
  }

  bool isprime()
    { return (roll[pos].isempty() || roll[pos].gettop()==n); }

  int64 nextprime()
  {
    do{
      next();
    }while(prevlen!=1);
    return prev[0];
  }

  ~Factgen() { delete [] roll; }


};

class FactgenWindow
{
    Factgen G;
    
public:
    int64 factors[60] ;
    
    int64* prev;
    int64* current;
    int64* next;
    
    int prevlen ;
    int64 prevval ;

    int currentlen ;
    int64 currentval ;
    
    int nextlen ;
    int64 nextval ;
    

    void init(int64 startin, int64 stopin)
    {
        G.init( startin, stopin);
        G.next();
        
        prev = &factors[0];
        current = &factors[20];
        next = &factors[40];
        
        memcpy(prev, G.prev, sizeof G.prev);
        
        prevlen = G.prevlen;
        prevval = G.prevn ;
        
        G.next();

        memcpy(current, G.prev, sizeof G.prev);
        currentlen = G.prevlen;
        currentval = G.prevn ;
        
        G.next();
        
        memcpy(next, G.prev, sizeof G.prev);
        nextlen = G.prevlen;
        nextval = G.prevn;
        
    }
    
    void advance()
    {
        G.next();
        
        int64* temp = prev;
        prev = current;
        current = next;
        next = temp;
        
        memcpy(next, G.prev, sizeof G.prev);
        
        prevlen = currentlen;
        prevval = currentval;
        
        currentlen = nextlen;
        currentval = nextval;

        nextlen = G.prevlen;
        nextval = G.prevn;
    }
    
    bool isprime()
    {
        return ( currentval == current[0] );
    }
      bool isprevprime()
    {
        return ( prevval == prev[0] );
    }
      bool isnextprime()
    {
        return ( nextval == next[0] );
    }
    
    bool issquarefree()
    {
        int64 prod = current[0];
        for( int i = 1; i < currentlen; i++)
        {
            prod = prod*current[i];
        }
        return ( prod == currentval );
    }
    
    bool isadmissible()
    {

      bool return_val = true;
      
      if(currentval%5==0)
      {
        return false;
      }
       if(currentval%2==0)
      {
        return false;
      }
      
      if(currentval == current[0])
      {
        return true;
      }
      else
      {
        int runningprod=current[0];

        for(int i=1; i<currentlen; i++)
        {
          if(currentval%5==1 || currentval%5==4)
          {
            if(gcd(currentval-1, runningprod)!=1)
            {
              return false;
            }

          }
          if(currentval%5==2 || currentval%5==3)
          {
            if(gcd(currentval+1, runningprod)!=1)
            {
              return false;
            }
          }
          runningprod*=current[i];

        }
        if( runningprod != currentval)
        {
          return false;
        }
      }
      return true;
    }

    
};



#endif
