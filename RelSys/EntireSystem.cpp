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
 * File:   EntireSystem.cpp
 * Author: Anders Reenberg Andersen
 * 
 */

#include "EntireSystem.h"
#include <iostream>
#include <math.h>

EntireSystem::EntireSystem(int nW, WardData* wards):
wards_pointer(wards),
nWards(nW)
{
    cout << "Initializing... " << flush;
    initializeSystem();
    cout << "Done." << endl;
}

EntireSystem::EntireSystem(const EntireSystem& orig) {
}

EntireSystem::~EntireSystem() {
}

void EntireSystem::initializeSystem(){
    
    //calculate local and global state space sizes
    nS = 1;
    for (int w=0; w<nWards; w++){
        calculateWardStateSpaceSize(w,nWards);
        nS *= getWardStateSpaceSize(w);
    }
    initializeState();
    initializeJumbVectors();
}

void EntireSystem::initializeState(){
    
    sidx = 0;
    state.clear(); state.resize(nWards); //state contains nWards different wards
    for (int i=0; i<nWards; i++){
        state[i].resize(nWards,0); //each ward contains nWards different patient types
    }
    K_use.clear(); K_use.resize(nWards,0);
    
}

void EntireSystem::initializeJumbVectors(){
    //initialize vectors that are used to track jumps (transitions)
    
    fromIdxSize = 0;
    
    //calculate maximum size
    int maximumSize = pow(nWards,2)*2;
    
    //allocate memory to vectors
    jumpFromIdx.clear(); jumpFromIdx.resize(maximumSize,-1);
    jumpFromRate.clear(); jumpFromRate.resize(maximumSize,-1);
    
}

void EntireSystem::nextCurrentState(){
    
    //advance state index
    if (sidx<(nS-1)){ 
        sidx++;
    }else{
        sidx = 0;
    }
    
    //advance state form
    int pidx,widx=0;
    do{
        pidx=nWards-1;
        do{
            if (K_use[widx]<getWardCapacity(widx)){
                state[widx][pidx]++; K_use[widx]++;
                pidx = -1; widx = nWards; //exit
            }else if (state[widx][pidx]>0){
                K_use[widx] -= state[widx][pidx];
                state[widx][pidx] = 0;
                pidx--;
            }else{
                pidx--;
            }
        }while (pidx>=0);
        widx++;
    }while (widx<nWards);
    
}

void EntireSystem::printState(){
    //prints the current state
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=0; pidx<nWards; pidx++){
            cout << state[widx][pidx] << " ";
        }
        cout << endl;
    }
    cout << endl; //make some space for the next print
    
}

void EntireSystem::allIngoing(){
    //generate all jump indices and rates relative to the current state
    
    fromIdxSize = 0;
    int jump, prod;
    vector<int> fwvec(nWards,0);
    
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=(nWards-1); pidx>=0; pidx--){
            
            prod = 1;
            
            //came from discharge
            if (K_use[widx]<getWardCapacity(widx)){
                //calculate jump
                if (widx>0){
                    for (int i=0; i<widx; i++){
                        prod *= getWardStateSpaceSize(i);
                    }
                }
                if (pidx<(nWards-1)){
                    jump = forwardOne(K_use[widx],fwvec,
                        (state[widx][pidx]+1),pidx,widx);
                }else{
                    jump = 1;
                }
                jump *= prod;
                jumpFromIdx[fromIdxSize] = sidx+jump;
                //calculate rate
                jumpFromRate[fromIdxSize] = getWardServiceRate(pidx)*(state[widx][pidx]+1);
                
                fromIdxSize++;
            }
            //came from admission
            if (state[widx][pidx]>0){
                //calculate jump
                if (prod==1 && widx>0){
                    for (int i=0; i<widx; i++){
                        prod *= getWardStateSpaceSize(i);
                    }
                }
                if (pidx<(nWards-1)){
                    jump = backwardOne(K_use[widx],fwvec,
                        (state[widx][pidx]-1),pidx,widx);
                }else{
                    jump = 1;
                }
                jump *= prod;
                jumpFromIdx[fromIdxSize] = sidx-jump;
                //calculate rate
                if (widx==pidx){
                    jumpFromRate[fromIdxSize] = getWardArrivalRate(pidx);
                }else{
                    jumpFromRate[fromIdxSize] = getWardArrivalRate(pidx)*
                            getWardRelocationProbabilities(pidx)[widx];
                }
                fromIdxSize++;
            }
        }
    }
    //calculate diagonal
    jumpFromIdx[fromIdxSize] = sidx;
    jumpFromRate[fromIdxSize] = diagonalRate();
    fromIdxSize++;
    
}

double EntireSystem::diagonalRate(){
    
    double diag = 0;
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=(nWards-1); pidx>=0; pidx--){
    
            //go to admission
            if (K_use[widx]<getWardCapacity(widx)){
                //calculate rate
                if (widx==pidx){
                    diag -= getWardArrivalRate(pidx);
                }else{
                    diag -= getWardArrivalRate(pidx)*
                            getWardRelocationProbabilities(pidx)[widx];
                }
            }
            //go to discharge
            if (state[widx][pidx]>0){
                //calculate rate
                diag -= getWardServiceRate(pidx)*state[widx][pidx];
            }
        }
    }

    return(diag);
}


int EntireSystem::forwardOne(int wardCapacityUsed, vector<int> &j, int targetval, int &pidx, int &widx){
    
    for (int i=0; i<nWards; i++){
        j[i] = state[widx][i];
    }
    
    int i, ii, m = 0;
    
    bool run = true;
    while(run){
        
        //advance state form by one step
        m++;
        
        i = nWards-1;
        while (i>=0){
            if (wardCapacityUsed<getWardCapacity(widx)){
                j[i]++; wardCapacityUsed++;
                i = -1; //exit
            }else if(j[i]>0){
                wardCapacityUsed -= j[i]; j[i] = 0;
                i--;
            }else{
                i--;
            }
        }
        
        if ( j[pidx]==targetval ){
            ii=0; run = false;
            do{   
                if (j[ii]!=state[widx][ii] && ii!=pidx){
                    run = true;
                }
                ii++;
            }while(run==false && ii<nWards);
            
        }
        
    }
    
    return(m);
}

int EntireSystem::backwardOne(int wardCapacityUsed, vector<int> &j, int targetval, int &pidx, int &widx){
    
    for (int i=0; i<nWards; i++){
        j[i] = state[widx][i];
    }
    j[pidx] = targetval;
    wardCapacityUsed--;
    
    int i, ii, m = 0;
    
    bool run = true;
    while(run){
        
        //advance state form by one step
        m++;
        
        i = nWards-1;
        while (i>=0){
            if (wardCapacityUsed<getWardCapacity(widx)){
                j[i]++; wardCapacityUsed++;
                i = -1; //exit
            }else if(j[i]>0){
                wardCapacityUsed -= j[i]; j[i] = 0;
                i--;
            }else{
                i--;
            }
        }
        
        if ( j[pidx]==state[widx][pidx] ){
            ii=0; run = false;
            do{   
                if (j[ii]!=state[widx][ii]){
                    run = true;
                }
                ii++;
            }while(run==false && ii<nWards);
            
        }
        
    }
    
    return (m);
}


int EntireSystem::getWardID(int ward){
    
    return((wards_pointer + ward)->wardnr);
}

double EntireSystem::getWardArrivalRate(int ward){
    
    return((wards_pointer + ward)->arrivalRate);
}

double EntireSystem::getWardServiceRate(int ward){
    
    return((wards_pointer + ward)->serviceRate);
}

vector<double> EntireSystem::getWardRelocationProbabilities(int ward){
    
    return((wards_pointer + ward)->relocationProbabilities);
}

void EntireSystem::calculateWardStateSpaceSize(int ward, int numberOfWards){
    
    return((wards_pointer + ward)->calculateWardStateSpace(numberOfWards));
}

int EntireSystem::getWardStateSpaceSize(int ward){
    
    return((wards_pointer + ward)->wardStateSpaceSize);
}

int EntireSystem::getWardCapacity(int ward){
    
    return((wards_pointer + ward)->capacity);
}