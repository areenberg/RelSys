/*
 * Copyright 2020 Anders Reenberg Andersen.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* 
 * File:   Combinatorial.cpp
 * Author: $Anders Reenberg Andersen
 * 
 * Created on 7. marts 2020, 23:41
 */

#include "Combinatorial.h"
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

Combinatorial::Combinatorial() {
    
}

Combinatorial::Combinatorial(const Combinatorial& orig) {
}

Combinatorial::~Combinatorial() {
}

//int combinatorial::capWeights(int K, int n) {
//    //derives the total number of permutations for a system
//    //with total capacity K and customer nodes with increasing weights.
//    //node weights start with 1 and increases with an
//    //increment of one. node state starts with zero customers and increases
//    //until the capacity is depleted.
//    
//    int ms = K+1;
//    if (n>1){
//        int m = 0, corr;
//        int j [(n-1)];
//        float lim;
//        int w [n]; //additional node weights
//        int ww = 2;
//        for (int i=0; i<(n-1); i++){
//            j[i] = 0;
//            w[i] = ww;
//            ww++;
//        }
//        
//        //run through all loops
//        bool run = true;
//        while(run){
//            
//            //adjust correction
//            corr = 0;
//            for (int i=0; i<(n-1); i++){ 
//                corr += w[i]*j[i];
//            }
//            m += (ms - corr); //increase sum
//        
//            //adjust j and limits
//            bool inc = true;
//            int inccount = 0;
//            for (int i=0; i<(n-1); i++){
//                //evaluate current limit of node
//                lim = K;
//                if (i<(n-2)){
//                    for (int ii=(i+1); ii<(n-1); ii++){
//                        lim -= w[ii]*j[ii];
//                    }
//                }
//                lim /= w[i];
//                lim = floor(lim);        
//                        
//                //adjust j
//                if (inc){
//                    if ( (j[i]+1) > lim ){
//                        j[i] = 0;
//                        inc = true;
//                        inccount++;
//                    }else{
//                        j[i]++;
//                        inc = false;
//                    }       
//                }
//            }
//            
//            //check if all loops are finished
//            if (inccount==(n-1)){
//                run = false;
//            }    
//        }
//        
//        return m;
//    }else{
//        return ms;
//    }
//}

int Combinatorial::capWithLimits(int &K, vector<int> &upperLim, vector<int> &lowerLim){
    //derives the total number of permutations for a system
    //with total capacity K and customer nodes corresponding to the
    //sizes of upperLim and lowerLim.
    //each node contributes to the capacity utilization with
    //a weight of one.
    //each node i can contain a minimum of lowerLim[i] items and a maximum
    //of upperLim[i] items.
    
    int j[upperLim.size()];
    K_use = 0;
    for (int i=0; i<upperLim.size(); i++){
        K_use += lowerLim[i];
        j[i] = lowerLim[i];
    }
    
    m = 1; //printJ(lng,j);
        
    int i = upperLim.size()-1;
    while (i>=0){
            
        if (i==(upperLim.size()-1)){
            while ( (j[i]+1)<=upperLim[i] && (K_use+1)<=K ){
                j[i]++; m++; K_use++; //printJ(lng,j);
            }
            K_use -= j[i]-lowerLim[i]; j[i] = lowerLim[i];
            i--;
        }else{
                
            if ( (j[i]+1)<=upperLim[i] && (K_use+1)<=K ){
                j[i]++; m++; K_use++; //printJ(lng,j);
                i = upperLim.size()-1;
            }else if(j[i]>lowerLim[i]){
                K_use -= j[i]-lowerLim[i]; j[i] = lowerLim[i];
                i--;
            }else{
                i--;
            }
                
        }
            
    }
        
    return (m);
}


//void combinatorial::printJ(int lng, int j[]){
//    
//    for (int i=0; i<lng; i++){
//        cout << j[i] << "  ";
//    }
//    cout << "" << endl;
//    
//}


//int combinatorial::capWithUpperBound(vector<int> &limits){
//    //corresponds to capWithoutWeights, but with an upper
//    //bound on the capacity in all nodes except for
//    //a "main node" which is able to contain the entire
//    //capacity. The methods finds this main node by looking
//    //for the node with the largest capacity.
//    
//    int n = limits.size();
//    int K = -1, idxmx = 0;
//    for (int i=0; i<n; i++){
//        if (limits[i]>K){
//            K = limits[i];
//            idxmx = i;
//        }
//    }
//    int lm[(n-1)]; //new limit array
//    int kk = 0;
//    for (int i=0; i<n; i++){
//        if (i!=idxmx){
//            lm[kk] = limits[i];
//            kk++;
//        }
//    }
//    
//    
//    int ms = K+1;
//    if (n>1){
//        int m = 0, corr;
//        int j [(n-1)];
//        float lim;
//        for (int i=0; i<(n-1); i++){
//            j[i] = 0;
//        }
//        
//        //run through all loops
//        bool run = true;
//        while(run){
//            
//            //adjust correction
//            corr = 0;
//            for (int i=0; i<(n-1); i++){ 
//                corr += j[i];
//            }
//            m += (ms - corr); //increase sum
//        
//            //adjust j and limits
//            bool inc = true;
//            int inccount = 0;
//            for (int i=0; i<(n-1); i++){
//                //evaluate current limit of node
//                lim = K;
//                if (i<(n-2)){
//                    for (int ii=(i+1); ii<(n-1); ii++){
//                        lim -= j[ii];
//                    }
//                }
//                
//                //user limits
//                if (lim>lm[i]){
//                    lim = lm[i];
//                }
//                
//                //adjust j
//                if (inc){
//                    if ( (j[i]+1) > lim ){
//                        j[i] = 0;
//                        inc = true;
//                        inccount++;
//                    }else{
//                        j[i]++;
//                        inc = false;
//                    }       
//                }
//            }
//            
//            //check if all loops are finished
//            if (inccount==(n-1)){
//                run = false;
//            }    
//        }
//        
//        return m;
//    }else{
//        return ms;
//    }
//    
//    return(0);
//}



//int combinatorial::capWithoutWeights_1(int K, int n){
//    //derives the total number of permutations for a system
//    //with total capacity K and customer nodes n.
//    //each node contributes to the capacity utilization with
//    //a weight of one.
//    //node state starts with zero customers and increases
//    //until the capacity is depleted.
//    
//    //-----------------------------------------------------
//    //SLOWER THAN capWithoutWeights_2, BUT STABLE FOR
//    //LARGE VALUES OF n.
//    //-----------------------------------------------------
//    
//    int ms = K+1;
//    if (n>1){
//        int m = 0, corr;
//        int j [(n-1)];
//        float lim;
//        for (int i=0; i<(n-1); i++){
//            j[i] = 0;
//        }
//        
//        //run through all loops
//        bool run = true;
//        while(run){
//            
//            //adjust correction
//            corr = 0;
//            for (int i=0; i<(n-1); i++){ 
//                corr += j[i];
//            }
//            m += (ms - corr); //increase sum
//        
//            //adjust j and limits
//            bool inc = true;
//            int inccount = 0;
//            for (int i=0; i<(n-1); i++){
//                //evaluate current limit of node
//                lim = K;
//                if (i<(n-2)){
//                    for (int ii=(i+1); ii<(n-1); ii++){
//                        lim -= j[ii];
//                    }
//                }
//                        
//                //adjust j
//                if (inc){
//                    if ( (j[i]+1) > lim ){
//                        j[i] = 0;
//                        inc = true;
//                        inccount++;
//                    }else{
//                        j[i]++;
//                        inc = false;
//                    }       
//                }
//            }
//            
//            //check if all loops are finished
//            if (inccount==(n-1)){
//                run = false;
//            }    
//        }
//        
//        return m;
//    }else{
//        return ms;
//    }
//    
//}


//int combinatorial::capWithoutWeights_2(int capacity, int nodes){
//    //derives the total number of permutations for a system
//    //with total capacity K and customer nodes.
//    //each node contributes to the capacity utilization with
//    //a weight of one.
//    //node state starts with zero customers and increases
//    //until the capacity is depleted.
//    
//    //-----------------------------------------------------
//    //FASTER THAN capWithoutWeights_1, BUT UNSTABLE FOR
//    //LARGE VALUES OF nodes.
//    //-----------------------------------------------------
//    
//    int m = 0;
//    for (int i=0; i<=capacity; i++){
//        m += numberAmongOwners(i,nodes);
//    }
//    
//    return m;
//}

//long double combinatorial::numberAmongOwners(int amount, int owners){
//    //calculates the number of permutations of distributing
//    //an amount of items among owners (for a single location).
//    
//    long double prod = 1;
//    int i;
//    for (i=1; i<=(owners-1); i++){
//        prod *= (amount + i);
//    }
//    
//    long double y = factorial((owners-1));
//    
//    prod *= (1/y);
//    
//    return prod;
//}


//long double combinatorial::factorial(int x){
//        //regular factorial
//    
//        long double fact = 1;
//         
//        if (x>0){
//            for (int i=1; i<=x; i++){
//                fact *= i;
//            }
//        }
//        
//        return fact;
//    }