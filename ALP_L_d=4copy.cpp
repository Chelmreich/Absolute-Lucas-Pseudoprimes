const int d = 7;
int jac_sym_lut[d] = {0, 1, 1, -1, 1, -1, -1};
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include "int.h"
#include "mpi.h"
#include "factgen.h"

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
    vector <int> primes_vec;
    vector <int> symbol_vec;
    primes_vec.push_back(3);
    symbol_vec.push_back(-1);
    //ofstream Large4("ALP_L_k=4_test.txt");
   
    string Name("L_k=4_");
    string rank = to_string(my_rank);
    string text(".txt");
    Name = Name + rank + text;
    ofstream Large4(Name); 

    //**** D=4 CASE ****
    int128 B=1;
    B = B<<64;
    int64 X=6000000;
    int64 i1=my_rank;
    int64 i2=0;
    int64 i3=0;
    int64 UPB1;
    int64 UPB2;
    int64 UPB3;
    UPB1=pow(B, (1.0/4));

    vector <int64> ALP_num;

    FactgenWindow P;
    int64 start = 5;
    int64 stop = 2000000;
    P.init(start,stop);

    while(P.currentval<=stop)
    {
        P.advance();

        if(P.isprime())
        {
            //put prime in vector
            primes_vec.push_back(P.currentval);
            symbol_vec.push_back(jac_sym_lut[P.currentval%d]);
        }
    }


    
    //auto begin = std::chrono::high_resolution_clock::now();
    Prime_init(primes_vec[i1], symbol_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);



    while(primes_vec[i1]<UPB1)
    {
        //cout << primes_vec[i1]  << endl;
        if(primes_vec[i1]>sqrt(X))
        {
            i2=i1+1;
        }
        else
        {
            i2=upper_bound(primes_vec.begin(), primes_vec.end(), (int) (X/Running_Prod[0]))-primes_vec.begin();
        }
        UPB2=pow((B/primes_vec[i1]), (1.0/3));
            
        //cout << i2  << " " << primes_vec[i2- 1]  << endl;


        while(primes_vec[i2]<UPB2)
        {
        
        
            if(gcd(primes_vec[i2]-symbol_vec[i2], Running_Prod[0])==1) //admissibility check
            {

                Prime_update(primes_vec[i2], symbol_vec[i2], 1, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 0);
                
                i3=i2+1;
                UPB3=pow((B/(primes_vec[i1]*primes_vec[i2])), (1.0/2));
                while(primes_vec[i3]<UPB3)
                {

                    if(gcd(primes_vec[i3]-symbol_vec[i3], Running_Prod[1])==1)
                    {
                        
                        Prime_update(primes_vec[i3], symbol_vec[i3], 2, Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign, 1);


                        int64 I=inv(Running_Prod[2], Running_L[2]);
                        int64 Pk_pos=I;
                        int64 Pk_neg = -I+Running_L[2];
                        
                        int128 PtimesL = (int128) Running_Prod[2] * (int128) Running_L[2];
                    
                        if(PtimesL>B)
                        {
                            if(Pk_pos>primes_vec[i3])
                            {
                                product = (int128) Running_Prod[2] * (int128) Pk_pos;
                                if(product<B)
                                {
                                    Primes_dividing_P[3] = Pk_pos;
                                    Sign[3] = jac_sym_lut[Pk_pos%d];
                                    if(KorseltCrit2(Primes_dividing_P, Sign, 4, product-jac_sym_lut[product%d])==true &&isPrime(Pk_pos))
                                    {
                                        Large4<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_pos<<endl;
                                        //cout<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_pos<<endl;
                                        //ALP_num.push_back(product);
                                    }
                                }
                            }

                            if(Pk_neg>primes_vec[i3])
                            {      
                                product = (int128) Running_Prod[2] * (int128)Pk_neg; 
                                if(product<B)
                                {
                                    Primes_dividing_P[3] = Pk_neg;
                                    Sign[3] = jac_sym_lut[Pk_neg%d];
                                    if(KorseltCrit2(Primes_dividing_P, Sign, 4, product-jac_sym_lut[product%d])==true &&isPrime(Pk_neg))
                                    {
                                        //ALP_num.push_back(product);
                                        Large4<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_neg<<endl;
                                        //cout<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_neg<<endl;
                                    }
                                }
                            }
                        }
                        else
                        {
                            
                            int64 K1=(Running_Prod[2]/Running_L[2])+1;

                            int64 K2=B/PtimesL+1;

                            int64 K=min(K1, K2);
                            for(int k=0; k<=K; k++)
                            {
                                
                        
                                if(Pk_pos>primes_vec[i3])
                                {
                                    product= (int128) Running_Prod[2] * (int128)Pk_pos;
                                                
                                    if(product<B)
                                    {
                                        Primes_dividing_P[3] = Pk_pos;
                                        Sign[3] = jac_sym_lut[Pk_pos%d];
                                        if(KorseltCrit2(Primes_dividing_P, Sign, 4, product-jac_sym_lut[product%d])==true &&isPrime(Pk_pos))
                                        {
                                            
                                            Large4<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_pos<<endl;
                                            //cout<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_pos<<endl;
                                            //ALP_num.push_back(product);
                                        }
                                    }
                                }

                                if(Pk_neg>primes_vec[i3])
                                {
                                    product= (int128) Running_Prod[2] * (int128)Pk_neg;
                                                    
                                    if(product<B)
                                    {
                                        Primes_dividing_P[3] = Pk_neg;
                                        Sign[3] = jac_sym_lut[Pk_neg%d];
                                        if(KorseltCrit2(Primes_dividing_P, Sign, 4, product-jac_sym_lut[product%d])==true &&isPrime(Pk_neg))
                                        {
                                            //ALP_num.push_back(product);
                                            Large4<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_neg<<endl;
                                            //cout<<product<<" "<<primes_vec[i1]<<" "<<primes_vec[i2]<<" "<<primes_vec[i3]<<" "<<Pk_neg<<endl;
                                        }
                                    }
                                }
                                Pk_pos = Pk_pos+Running_L[2];
                                Pk_neg = Pk_neg+Running_L[2];
                            }
                        }
                    }
                    i3++;
                }
            }
            i2++;
        }
        i1+=proc;
        //i1++;
        Prime_init(primes_vec[i1], symbol_vec[i1], Primes_dividing_P, Prime_symbol, Sign, Running_L, Running_Prod, Running_Sign);
        
    }
    
    MPI_Finalize();

    /*auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    

    
    printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    cout<<ALP_num.size()<<endl;*/


    Large4.close();
    
        
}
    
