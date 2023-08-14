#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include "int.h"
//#include "mpi.h"
#include "factgen.h"

using namespace std;

int main(int argc, char * argv[])
{
    /*int my_rank;
    int proc;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &proc);*/

    int64 Primes_dividing_P[12] = {};
    int64 Prime_symbol[12] = {};
    int64 Sign[12] = {};
    int64 Running_L[12] = {};
    int64 Running_Prod[12] = {};
    int64 Running_Sign[12] = {};
    int128 product;
    vector <int> primes_vec;
    primes_vec.push_back(3);

    //ofstream Large5("Large_k=5.txt");

    /*string Name("L_k=5_");
    string rank = to_string(my_rank);
    string text(".txt");
    Name = Name + rank + text;
    ofstream Large5(Name);*/

    int128 B=10000000000;
    //B = B<<10;
    //int64 X=6000000;
    int64 i1=0;
    int64 i2=0;
    int64 i3=0;
    int64 i4=0;
    int64 UPB1;
    int64 UPB2;
    int64 UPB3;
    int64 UPB4;
    UPB1=pow(B, (1.0/5));

    vector <int64> ALP_num;

    FactgenWindow P;
    int64 start = 5;
    int64 stop = 2000000;
    P.init(start, stop);

    while(P.currentval<=stop)
    {
        P.advance();

        if(P.isprime())
        {
            primes_vec.push_back(P.currentval);
        }
    }


    auto begin = std::chrono::high_resolution_clock::now();
    Prime_init(primes_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);

    while(primes_vec[i1]<UPB1)
    {
        i2=i1+1;
        
        UPB2=pow((B/Running_Prod[0]), (1.0/4));

        while(primes_vec[i2]<UPB2)
        {
            if(gcd(primes_vec[i2]-jacobi(5, primes_vec[i2]), Running_Prod[0])==1)
            {
                Prime_update(primes_vec[i2], 1, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 0);

                /*if(primes_vec[i2]>sqrt(X/Running_Prod[0]))
                {
                    i3=i2+1;
                }
                else
                {
                    i3=upper_bound(primes_vec.begin(), primes_vec.end(), (int) (X/Running_Prod[1]))-primes_vec.begin();
                }*/
                i3=i2+1;
                UPB3=pow((B/Running_Prod[1]), (1.0/3));

                while(primes_vec[i3]<UPB3)
                {
                    if(gcd(primes_vec[i3]-jacobi(5, primes_vec[i3]), Running_Prod[1])==1)
                    {
                        Prime_update(primes_vec[i3], 2, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 1);

                        i4=i3+1;
                        UPB4=pow((B/Running_Prod[2]), (1.0/2));

                        while(primes_vec[i4]<UPB4)
                        {
                            if(gcd(primes_vec[i4]-jacobi(5, primes_vec[i4]), Running_Prod[2])==1)
                            {
                                Prime_update(primes_vec[i4], 3, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 2);

                                int64 I=inv(Running_Prod[3], Running_L[3]);

                                int64 K1=(Running_Prod[3]/Running_L[3])+1;

                                int64 K2=(B/Running_Prod[3]*Running_L[3])+1;

                                int64 K=min(K1, K2);

                                for(int k=0; k<K; k++)
                                {
                                    int64 Pk_pos=I+(k*Running_L[3]);
                                    int64 Pk_neg=-I+(k*Running_L[3]);

                                    if(Pk_pos>primes_vec[i1] && Pk_pos>primes_vec[i2] && Pk_pos>primes_vec[i3] && Pk_pos>primes_vec[i4])
                                    {
                                        if(KorseltCrit(Primes_dividing_P, Running_Prod[2], 3, primes_vec[i4], jacobi(5, primes_vec[i4]), jacobi(5, Pk_pos), Pk_pos, jacobi(5, Running_Prod[2]))==true && isPrime(primes_vec[i2]) && isPrime(primes_vec[i3]) && isPrime(primes_vec[i4]) &&isPrime(Pk_pos))
                                        {
                                            product= (int128) primes_vec[i1] * (int128)primes_vec[i2]* (int128)primes_vec[i3]* (int128)primes_vec[i4]* (int128)Pk_pos;
                                            //Large5<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_pos<<endl;
                                            ALP_num.push_back(product);
                                        }
                                    }

                                    if(Pk_neg>primes_vec[i1] && Pk_neg>primes_vec[i2] && Pk_neg>primes_vec[i3] && Pk_neg>primes_vec[i4])
                                    {
                                        if(KorseltCrit(Primes_dividing_P, Running_Prod[2], 3, primes_vec[i4], jacobi(5, primes_vec[i4]), jacobi(5, Pk_neg), Pk_neg, jacobi(5, Running_Prod[2]))==true && isPrime(primes_vec[i2]) && isPrime(primes_vec[i3]) && isPrime(primes_vec[i4]) &&isPrime(Pk_neg))
                                        {
                                            product=(int128)primes_vec[i1]*(int128)primes_vec[i2]*(int128)primes_vec[i3]*(int128)primes_vec[i4]*(int128)Pk_neg;
                                            //Large5<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_neg<<endl;
                                            ALP_num.push_back(product);
                                        }
                                    }
                                   

                                        cout<<"Pk: "<<Pk_neg<<"L: "<<Running_L[3]<<" K: "<<k<<endl;
                                        cout<<"Pk: "<<Pk_pos<<"L: "<<Running_L[3]<<" K: "<<k<<endl;

                                }
                            } 
                            i4++;
                            
                        }
                
                    }
                    i3++;
                    
                }
                
            }
            i2++;
        }
        //i1+=proc;
        i1++;
        Prime_init(primes_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);
    }
    //MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    //printf("Result: %.20f\n", sum);
    
    //printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);

    //cout<<ALP_num.size()<<endl;
    Large5.close();
}