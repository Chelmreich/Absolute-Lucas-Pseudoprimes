const int d = 7;
int jac_sym_lut[d] = {0, 1, 1, -1, 1, -1, -1};
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include "int.h"
#include "mpi.h"
#include "factgen.h"
// mpicxx -o blah file.cpp
// mpirun -np 32 blah
using namespace std;

int main(int argc, char * argv[])
{
    int my_rank;
    int proc;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    int64 Primes_dividing_P[12] = {};
    int64 Prime_symbol[12] = {};
    int64 Sign[12] = {};
    int64 Running_L[12] = {};
    int64 Running_Prod[12] = {};
    int64 Running_Sign[12] = {};
    int128 product;
    int128 product_neg;
    int128 product_pos;
    vector <int> primes_vec;
    vector <int> symbol_vec;
    //vector <int> primeminussymbol_vec;
    primes_vec.push_back(3);
    symbol_vec.push_back(-1);
    

    //ofstream Large5("Large_k=5_compare.txt");

    string Name("L_k=5_");
    string rank = to_string(my_rank);
    string text(".txt");
    Name = Name + rank + text;
    ofstream Large5(Name);

    int128 B=1;
    B = B<<64;
    int64 X=6000000;
    int64 i1=my_rank;
    int64 i2=0;
    int64 i3=0;
    int64 i4=0;
    int64 UPB1;
    int64 UPB2;
    int64 UPB3;
    int64 UPB4;
    int64 work=0;

    

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
            symbol_vec.push_back(jac_sym_lut[P.currentval%d]);
            //primeminussymbol_vec.push_back(P.currentval-jac_sym_lut[P.currentval%5]);
        }
    }



    //auto begin = std::chrono::high_resolution_clock::now();
    Prime_init(primes_vec[i1], symbol_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);
    
    UPB1=pow(B, (1.0/5));
    
    while(primes_vec[i1]<UPB1)
    {
        //cout<<primes_vec[i1]<<endl;
        //Prime_init(primes_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);
        i2=i1+1;
        
        UPB2=pow((B/Running_Prod[0]), (1.0/4));

       
        while(primes_vec[i2]<UPB2)
        {
            //cout<<primes_vec[i2]<<endl;
            if(gcd(primes_vec[i2]-symbol_vec[i2], Running_Prod[0])==1)
            {
                Prime_update(primes_vec[i2], symbol_vec[i2], 1, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 0);
                if(primes_vec[i1]>sqrt(X/Running_Prod[0]))
                {
                    i3=i2+1;
                }
                else
                {
                    i3=upper_bound(primes_vec.begin(), primes_vec.end(), (int) (X/Running_Prod[1]))-primes_vec.begin();
                }
                //i3=i2+1;
                UPB3=pow((B/Running_Prod[1]), (1.0/3));

                work++;
                work = work % proc;
                while(my_rank == work && primes_vec[i3]<UPB3)
                //while(primes_vec[i3]<UPB3)
                {
                    //cout<<primes_vec[i3]<<endl;
                    if(gcd(primes_vec[i3]-symbol_vec[i3], Running_Prod[1])==1)
                    {
                        Prime_update(primes_vec[i3], symbol_vec[i3], 2, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 1);

                        i4=i3+1;
                        UPB4=pow((B/Running_Prod[2]), (1.0/2));
                        //cout<<UPB4<<" "<<Running_Prod[2]<<endl;
                        while(primes_vec[i4]<UPB4)
                        {
                            //cout<<primes_vec[i4]<<endl;
                            if(gcd(primes_vec[i4]-symbol_vec[i4], Running_Prod[2])==1)
                            {
                                Prime_update(primes_vec[i4], symbol_vec[i4], 3, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 2);

                                int64 I=inv(Running_Prod[3], Running_L[3]);
                                int64 Pk_pos=I;
                                int64 Pk_neg = -I+Running_L[3];
                                        
                                int128 PtimesL = (int128) Running_Prod[3] * (int128) Running_L[3];
                                if(PtimesL>B)
                                {
                                    if(Pk_pos>primes_vec[i4])
                                    {
                                        product_pos=(int128)Running_Prod[3]*(int128)Pk_pos;
                                        if(product_pos<B)
                                        {
                                            Primes_dividing_P[4] = Pk_pos;
                                            Sign[4] = jac_sym_lut[Pk_pos%d];
                                            if(KorseltCrit2(Primes_dividing_P, Sign, 5, product_pos-jac_sym_lut[product_pos%d])==true &&isPrime(Pk_pos))
                                            {
                                                Large5<<product_pos<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_pos<<endl;
                                                //ALP_num.push_back(product_pos);
                                            }
                                        }
                                    }

                                    if(Pk_neg>primes_vec[i4])
                                    {
                                        product_neg=(int128)Running_Prod[3]*(int128)Pk_neg;
                                        if(product_neg<B)
                                        {
                                            Primes_dividing_P[4] = Pk_neg;
                                            Sign[4] = jac_sym_lut[Pk_neg%d];
                                            if(KorseltCrit2(Primes_dividing_P, Sign, 5, product_neg-jac_sym_lut[product_neg%d])==true &&isPrime(Pk_neg))
                                            {
                                                Large5<<product_neg<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_neg<<endl;
                                                //ALP_num.push_back(product_neg);
                                            }
                                        }
                                    }    
                                }
                                
                                else
                                {
                                    
                                    int64 K1=(Running_Prod[3]/Running_L[3])+1;

                                    int64 K2 = B/PtimesL + 1;

                                    int64 K=min(K1, K2);
                                    for(int k=0; k<=K; k++)
                                    {
                                        if(Pk_pos>primes_vec[i4])
                                        {
                                            product= (int128) Running_Prod[3]* (int128)Pk_pos;
                                            if(product<B)
                                            {
                                                Primes_dividing_P[4] = Pk_pos;
                                                Sign[4] = jac_sym_lut[Pk_pos%d];
                                                if(KorseltCrit2(Primes_dividing_P, Sign, 5, product-jac_sym_lut[product%d])==true &&isPrime(Pk_pos))
                                                {
                                                    Large5<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_pos<<endl;
                                                    //ALP_num.push_back(product);
                                                }
                                            }
                                        }

                                        if(Pk_neg>primes_vec[i4])
                                        {
                                            product=(int128)Running_Prod[3]*(int128)Pk_neg;
                                            if(product<B)
                                            {
                                                Primes_dividing_P[4] = Pk_neg;
                                                Sign[4] = jac_sym_lut[Pk_neg%d];
                                                if(KorseltCrit2(Primes_dividing_P, Sign, 5, product-jac_sym_lut[product%d])==true &&isPrime(Pk_neg))
                                                {
                                                    Large5<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<primes_vec[i4]<<" "<<Pk_neg<<endl;
                                                    //ALP_num.push_back(product);
                                                }
                                            }
                                        }
                                        //cout<<"Pk: "<<Pk_neg<<"L: "<<Running_L[3]<<" K: "<<k<<endl;
                                        //cout<<"ALP: "<<product<<" Pk: "<<Pk_pos<<" L: "<<Running_L[3]<<" K: "<<k<<endl;
                                        Pk_pos=Pk_pos+Running_L[3];
                                        Pk_neg=Pk_neg+Running_L[3];
                                    }
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
        Prime_init(primes_vec[i1], symbol_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);
    }
    MPI_Finalize();

    /*auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    

    
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);*/

    //sort(ALP_num.begin(), ALP_num.end());
    /*for(int i=0; i<ALP_num.size(); i++)
    {
        cout<<ALP_num[i]<<endl;
    }*/
    //cout<<ALP_num.size()<<endl;
    /*for(int i=0; i<ALP_num.size(); i++)
    {
       Large5<<ALP_num[i]<<endl; 
    }*/
    
    //cout<<i1<<endl;
    Large5.close();
}
