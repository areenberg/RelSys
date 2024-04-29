
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

int Combinatorial::capWithLimits(int &K, vector<int> &upperLim, vector<int> &lowerLim){
    //derives the total number of permutations for a system
    //with total capacity K and customer nodes corresponding to the
    //sizes of upperLim and lowerLim.
    //each node contributes to the capacity utilization with
    //a weight of one.
    //each node i can contain a minimum of lowerLim[i] items and a maximum
    //of upperLim[i] items.
    
    vector<int> j(upperLim.size());
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
