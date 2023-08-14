const int d =5;
int jac_sym_lut[d] = {0,1,-1,-1,1};

#include <iostream>
#include <fstream>
#include <vector>
#include "int.h"
#include "factgen.h"
#include "mpi.h"
#include <string.h>
#include <chrono>
using namespace std;


int main(int argc, char * argv[])
{
    int64 ub;
    int64 Exp_dP[20] = {};
    int64 Exp_PD[20] = {};
    int64 prim_merge[20] = {};
    int64 exp_merge[20] = {};


    int64 length;

    vector <int64> divisor;
    int64 divlen;
    int64 D;
    int64 C;
    int128 B=1;
    B = B<<64;
    int64 r;
    int64 q;
    int64 lrg_prim_fac;
    int64 eq;
    int64 er;
    int64 eP;
    int128 n; 
    //ofstream ALP("timingD3.txt");

    

    vector <int64> ALP_num;

    int kvalues[10] = {};
    int64 *deltaP;
    int64 delta_len;
    int64 delta_val;
    int64 count=0;
    int my_rank;
    int proc;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    string Name3("k=3_");
    string rank3 = to_string(my_rank);
    string text3(".txt");
    Name3 = Name3 + rank3 + text3;
    ofstream output3(Name3);

    string Name4("k=4_");
    string rank4 = to_string(my_rank);
    string text4(".txt");
    Name4 = Name4 + rank4 + text4;
    ofstream output4(Name4);    

    string Name5("k=5_");
    string rank5 = to_string(my_rank);
    string text5(".txt");
    Name5 = Name5 + rank5 + text5;
    ofstream output5(Name5);

    string Name6("k=6_");
    string rank6 = to_string(my_rank);
    string text6(".txt");
    Name6 = Name6 + rank6 + text6;
    ofstream output6(Name6);

    string Name7("k=7_");
    string rank7 = to_string(my_rank);
    string text7(".txt");
    Name7 = Name7 + rank7 + text7;
    ofstream output7(Name7);
    string Name("CD_off_4");
    string rank = to_string(my_rank);
    string text(".txt");
    Name = Name + rank + text;
    ofstream ALP(Name);

    FactgenWindow P;
    int64 start = 1;
    int64 stop =  6000000;
    P.init(start,stop);

    int64 work = 0;
    //auto begin = std::chrono::high_resolution_clock::now();

    while(P.currentval<=stop)
    {
        P.advance();

        if( P.isboundsadmissible(B) && P.isadmissible() )
        {
            
            //cout<<P.currentval<<endl;

            work++;
            if( (work % proc) == my_rank )
            {

              FactgenWindow PplusD;
              ub = 2*P.currentval-1;
              PplusD.init(1, ub);

              if(jac_sym_lut[ P.currentval % d] == 1 )
              {
                deltaP = P.prev;
                delta_len = P.prevlen;
                delta_val = P.prevval;
                eP = 1;
              }
              else
              {
                deltaP = P.next;
                delta_len = P.nextlen;
                delta_val = P.nextval;
                eP = -1;
              }

              // create exponent array for P - d_P
              createExpArr( Exp_dP, deltaP, delta_len, delta_val);
              // can decrement the exponent on 2 to guarantee q will be odd
              Exp_dP[0]--;
              // storing the exponent so we can restore it
              int64 exp_on_2 = Exp_dP[0];
              // we'll need the largest prime factor of P
              lrg_prim_fac = P.current[ P.currentlen - 1 ];

              while( PplusD.prevval<=ub )
              {
                
                // D up to sign
                D=PplusD.prevval-P.currentval;
                // fix the sign
                if( D < 0 ) { D = -D; }
                // D = 0 is not allowed
                if( D == 0 ){ PplusD.advance(); continue;}

                // need C to be integral, D is even, means divsor needs to be odd, too
                if( D % 2 == 0 )
                {
                    // no powers of 2
                    Exp_dP[0] = 0;
                }
                else
                {
                  // restore to default
                   Exp_dP[0] = exp_on_2;
                }
                //create exponent array of PplusD
                createExpArr( Exp_PD, PplusD.prev, PplusD.prevlen, PplusD.prevval);
                //merges the prime factorizations
                Merge(deltaP, PplusD.prev, Exp_dP, Exp_PD, prim_merge, exp_merge, delta_len, PplusD.prevlen, length);
                //create the divisors
                getDivisorsVec(prim_merge, exp_merge, 1, length, divisor);

                //cout<< P.currentval << " " << divisor.capacity() << " " << divisor.size() << " " << PplusD.prevval << " " << D << endl;

                //eq and er have opposite signs
                if(PplusD.prevval<P.currentval)
                {
                  for( int i = 0; i < divisor.size(); i++)
                  {
                    for( int eq = -1; eq < 2; eq +=2 )
                    {
                      er = -eq;
                      if( ( ( P.currentval*P.currentval ) + eq*divisor[i] ) % D == 0 )
                      {
                          C = ( ( ( P.currentval*P.currentval ) + eq*divisor[i] ) / D );
                          q = ( ( delta_val * PplusD.prevval ) / divisor[i] ) + eq;
                          if( jac_sym_lut[ q % d ] == eq
                            && q > lrg_prim_fac
                            && ( (int128)delta_val * (int128)( -P.currentval + C ) ) % (int128)divisor[i]==0
                            && ((int128)P.currentval*(int128)q*(int128)q) < B )
                          {
                              r =(((int128)delta_val*(int128)(-P.currentval + C))/(int128)divisor[i]) + er;
                              if( jac_sym_lut[r%d] == er
                                  && ((int128)P.currentval*(int128)q*(int128)r) < B
                                  && KorseltCrit(P.current, P.currentval, P.currentlen, q, eq, er, r, eP)
                                  && isPrime(r)
                                  && isPrime(q) )
                              {
                                  n = (int128)P.currentval*(int128)q*(int128)r;
                                  ALP<<n<<endl;
                                  switch(P.currentlen){
                                    case 1 : 
                                        output3 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output3<< P.current[i] << " " ;
                                        }
                                        output3 << q << " " << r << " " << endl;
                                        break;
                                    case 2 : 
                                        output4 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output4<< P.current[i] << " " ;
                                        }
                                        output4 << q << " " << r << " " << endl;
                                        break;
                                    case 3 : 
                                        output5 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output5<< P.current[i] << " " ;
                                        }
                                        output5 << q << " " << r << " " << endl;
                                        break;
                                    case 4 : 
                                        output6 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output6<< P.current[i] << " " ;
                                        }
                                        output6 << q << " " << r << " " << endl;
                                        break;
                                    case 5 : 
                                        output7 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output7<< P.current[i] << " " ;
                                        }
                                        output7 << q << " " << r << " " << endl;
                                        break;
                                    default :
                                        cout<<"oops"<<endl;
                                        cout<< P.currentval << " " << q << " " << r << endl;

                                  }

                            


                                  //ALP<<P.currentval<<" "<<q<<" "<<r<<endl;
                                  //ALP_num.push_back(P.currentval);
                                  //cout<<P.currentval<<" "<<q<<" "<<r<<endl;
                              }
                            }
                        }
                      }
                    }
                  }
                // eq and er have the same signs
                if(PplusD.prevval>P.currentval)
                {
                  for(int i=0; i<divisor.size(); i++)
                  {
                    for( int eq = -1; eq < 2; eq +=2 )
                    {
                      er = eq;
                      if( ( ( P.currentval*P.currentval ) + eq*divisor[i] ) % D == 0 )
                      {
                        C = ( ( ( P.currentval*P.currentval ) + eq*divisor[i] ) / D );
                        q = ( ( delta_val * PplusD.prevval ) / divisor[i] ) + eq;

                        if( jac_sym_lut[ q % d ] == eq
                          && q > lrg_prim_fac
                          && ( (int128)delta_val * (int128)( P.currentval + C ) ) % (int128)divisor[i] == 0 
                          && ((int128)P.currentval*(int128)q*(int128)q) < B )
                          {
                            r =( ( (int128)delta_val* (int128)( P.currentval + C ) )/ (int128)divisor[i] ) + er;
                            if( jac_sym_lut[r%d] == er
                                && ((int128)P.currentval*(int128)q*(int128)r) < B
                                && KorseltCrit(P.current, P.currentval, P.currentlen, q, eq, er, r, eP)
                                && isPrime(r)
                                && isPrime(q) )
                            {

                                n = (int128)P.currentval*(int128)q*(int128)r;
                                ALP<<n<<endl;
                                  switch(P.currentlen){
                                    case 1 : 
                                        output3 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output3<< P.current[i] << " " ;
                                        }
                                        output3 << q << " " << r << " " << endl;
                                        break;
                                    case 2 : 
                                        output4 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output4<< P.current[i] << " " ;
                                        }
                                        output4 << q << " " << r << " " << endl;
                                        break;
                                    case 3 : 
                                        output5 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output5<< P.current[i] << " " ;
                                        }
                                        output5 << q << " " << r << " " << endl;
                                        break;
                                    case 4 : 
                                        output6 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output6<< P.current[i] << " " ;
                                        }
                                        output6 << q << " " << r << " " << endl;
                                        break;
                                    case 5 : 
                                        output7 << n << " " ;
                                        for( int i=0; i < P.currentlen; i++)
                                        {
                                            output7<< P.current[i] << " " ;
                                        }
                                        output7 << q << " " << r << " " << endl;
                                        break;
                                    default :
                                        cout<<"oops"<<endl;
                                        cout<< P.currentval << " " << q << " " << r << endl;

                                  }
                                //ALP<<P.currentval<<" "<<q<<" "<<r<<endl;
                                //ALP_num.push_back(P.currentval);

                                //cout<<P.currentval<<" "<<q<<" "<<r<<endl;
                            }
                          }
                        }
                      }
                    }
                  }
                PplusD.advance();
              }
            }
        }
    }
    /*auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);*/
    //cout<<ALP_num.size()<<endl;

    MPI_Finalize();
    output3.close();
    output4.close();
    output5.close();
    output6.close();
    output7.close();
    ALP.close();
}
