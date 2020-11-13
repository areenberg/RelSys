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
#include "StatusBar.h"
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


void EntireSystem::nextStateLocalInput(vector<vector<int>> &s, vector<int> &capuse){
    //advance local input state form
    
    int pidx,widx=0;
    do{
        pidx=nWards-1;
        do{
            if (capuse[widx]<getWardCapacity(widx)){
                s[widx][pidx]++; capuse[widx]++;
                pidx = -1; widx = nWards; //exit
            }else if (s[widx][pidx]>0){
                capuse[widx] -= s[widx][pidx];
                s[widx][pidx] = 0;
                pidx--;
            }else{
                pidx--;
            }
        }while (pidx>=0);
        widx++;
    }while (widx<nWards);
    
}

void EntireSystem::selfCheck(){
    cout << "\nRUNNING SELF-CHECK" << endl;
    
    //store state space
    cout << "Storing state space..." << flush;
    vector<vector<vector<int>>> stateSpace(nS);
    initializeState();
    for (int i=0; i<nS; i++){
        
        stateSpace[i].resize(nWards);
        for (int j=0; j<nWards; j++){
            stateSpace[i][j].resize(nWards,0);
        }
        for (int j=0; j<nWards; j++){
            for (int k=0; k<nWards; k++){
                stateSpace[i][j][k] = state[j][k];
            }
        }
        nextCurrentState();
    }
    cout << " done." << endl;
    
    bool allOk; //indicates if errors were found
    
    cout << "Checking on-demand column indices and rates:" << endl;
    initializeState();
    for (int i=0; i<nS; i++){
        allIngoing();
        
        for (int j=0; j<fromIdxSize; j++){
            int idx = jumpFromIdx[j];
            double rate = jumpFromRate[j];
            allOk = checkOneStateChange(i,stateSpace[i],stateSpace[idx],rate,idx);
        }
        
        nextCurrentState();
    }
    cout << "Done." << endl;
    if (allOk){
        cout << "No errors were detected." << endl;
    }
    
    cout << "Checking stored transition matrix:" << endl;
    cout << "Building transposed transition matrix." << endl;
    buildTransposedChain();
    for (int i=0; i<nS; i++){
        if (qColumnIndices[i].size()!=qValues[i].size()){
            cout << "qColumnIndices["<<i<<"] has different size than qValues["<<i<<"]"<<endl;
            allOk = false;
        }
        for (int j=0; j<qColumnIndices[i].size(); j++){
            int idx = qColumnIndices[i][j];
            double rate = qValues[i][j];
            allOk = checkOneStateChange(i,stateSpace[i],stateSpace[idx],rate,idx);
        }
    }
    //check diagonal
    cout << "Back-transposing to actual transition matrix..." << flush;
    transposeChain();
    cout << " done." << endl;
    allOk = checkDiagonalRates();
    cout << "Done." << endl;
    if (allOk){
        cout << "No errors were detected." << endl;
    }
    
    cout << "SELF-CHECK COMPLETED \n" << endl;
    
}

bool EntireSystem::checkOneStateChange(int currentStateIdx, vector<vector<int>> &sc, vector<vector<int>> &sf, double rate, int idx){
    //self-check method
    
    bool allOk = true;
    
    //state differences
    int widx, pidx, dis = 0, adm = 0;
    double diff, diff1, targetRate;
    for (int i=0; i<sf.size(); i++){
        for (int j=0; j<sf[i].size(); j++){
            diff = sf[i][j]-sc[i][j];
            if (diff>0){ //discharge
                dis++;
                widx = i; pidx = j;
                diff1 = diff;
            }else if(diff<0){
                adm++;
                widx = i; pidx = j;
                diff1 = diff;
            }
        }
    }
    
    if ( (adm+dis)>1 ){
        cout << "ERROR: Wrong state change. Discharges: " << dis << ", Admissions: " << adm << endl;
        cout << "Current state index: " << currentStateIdx << ", jump from index: " << idx << endl;
        allOk = false;
    }
    
    //check rate
    if ( (adm+dis)==1 ){
        if (diff1<0){
            if (widx==pidx){
                targetRate = getWardArrivalRate(widx);
            }else{
                targetRate = getWardArrivalRate(pidx)*getWardRelocationProbabilities(pidx)[widx];
            }
            if (targetRate!=rate){
                cout << "ERROR: Wrong arrival rate. Target rate: " << targetRate << ", actual rate: " << rate << endl;
                cout << "Current state index: " << currentStateIdx << ", jump from index: " << idx << endl;
                cout << "Ward index: " << widx << ", patient index: " << pidx << endl;
                allOk = false;
            }
        }else{
            targetRate = getWardServiceRate(pidx)*sf[widx][pidx];
            if (targetRate!=rate){
                cout << "ERROR: Wrong service rate. Target rate: " << targetRate << ", actual rate: " << rate << endl;
                cout << "Current state index: " << currentStateIdx << ", jump from index: " << idx << endl;
                cout << "Ward index: " << widx << ", patient index: " << pidx << endl;
                cout << "Patients in from state: " << sf[widx][pidx] << ", patients in current state: " << sc[widx][pidx] << endl;
                allOk = false;
            }
        }
    }else if( (adm+dis)==0 ){
        if (rate>=0){
            cout << "ERROR: Diagonal rate appears to be non-negative." << endl;
            cout << "adm+dis = " << (adm+dis) << ", rate = " << rate << endl;
            allOk = false;
        }
    }
    
    return(allOk);
    
}

void EntireSystem::transposeChain(){
    //self-check method
    //transpose the stored transition matrix and
    //store the result in qColumnIndices_t and qValues_t.
    
    //allocate memory
    qColumnIndices_t.resize(qColumnIndices.size());
    qValues_t.resize(qValues.size());
    
    int m,j;
    
    //build
    for (int s=0; s<qColumnIndices.size(); s++){
        //count
        m = 1;
        for (int i=0; i<qColumnIndices.size(); i++){
            j = 0;
            while(s!=qColumnIndices[i][j] && j<(qColumnIndices[i].size()-1)){
                j++;
            }
            if (s==qColumnIndices[i][j] && j<(qColumnIndices[i].size()-1)){
                m++;
            }
        }
        
        //write values
        qColumnIndices_t[s].resize(m,0);
        qValues_t[s].resize(m,0);
        qColumnIndices_t[s][m-1] = s;
        qValues_t[s][m-1] = qValues[s][qValues[s].size()-1];
        m = 0;
        for (int i=0; i<qColumnIndices.size(); i++){
            j = 0;
            while(s!=qColumnIndices[i][j] && j<(qColumnIndices[i].size()-1)){
                j++;
            }
            if (s==qColumnIndices[i][j] && j<(qColumnIndices[i].size()-1)){
                qColumnIndices_t[s][m] = i;
                qValues_t[s][m] = qValues[i][j];
                m++;
            }
        }
        
    }
    
}

bool EntireSystem::checkDiagonalRates(){
    //self-check method.
    //checks if the value of the diagonal rate
    //equals the negative sum of the off-diagonal rates
    //in each row of the transition matrix.
    
    bool allOk = true;
    double sm, diag;
    for (int s=0; s<nS; s++){
        sm = 0;
        for (int j=0; j<(qValues_t[s].size()-1); j++){
            sm += qValues_t[s][j];
        }
        sm *= -1;
        diag = qValues_t[s][qValues_t[s].size()-1];
        double diff = sm-diag;
        if (abs(diff)>1e-9){
            allOk = false;
            cout << "ERROR: Diagonal rate does not equal the negative sum of the off-diagonal rates." << endl;
            cout << "State index: " << s << ", diagonal rate " << diag << ", should have equaled: " << sm << ". Difference: " << diff << endl;    
        }
    }
    
    return(allOk);
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

void EntireSystem::buildTransposedChain(){
    //creates and stores the entire <transposed> transition rate matrix
    
    qColumnIndices.clear(); qColumnIndices.resize(nS);
    qValues.clear(); qValues.resize(nS);
    initializeState();
    
    StatusBar sbar(nS,30);
    for (int i=0; i<nS; i++){
        allIngoing();
        qColumnIndices[i].resize(fromIdxSize,0);
        qValues[i].resize(fromIdxSize,0);
        for (int j=0; j<fromIdxSize; j++){
            qColumnIndices[i][j] = jumpFromIdx[j];
            qValues[i][j] = jumpFromRate[j];   
        }
        nextCurrentState();
        if (i%100000==0){
            double bval = i;
            sbar.updateBar(bval);
        }
    }
    sbar.endBar();
    
}

vector<double> EntireSystem::marginalDist(int widx, vector<double> &pi){
    
    vector<double> dist((getWardCapacity(widx)+1),0);
    
    int occupancy;
    initializeState();
    for (int i=0; i<nS; i++){
        occupancy = 0;
        for (int j=0; j<state[widx].size(); j++){
            occupancy += state[widx][j];
        }
        dist[occupancy] += pi[i];
        nextCurrentState();
    }
    
    return(dist);
    
}

double EntireSystem::expectedOccupancy(int widx, vector<double> &pi){
    
    vector<double> margDist = marginalDist(widx,pi);
    
    double y = 0; int k = 0;
    for (int i=0; i<=getWardCapacity(widx); i++){
        y += i*margDist[k];
        k++;
    }
    
    return(y);
}

double EntireSystem::expectedOccupancyFraction(int widx, vector<double> &pi){
    
    return(expectedOccupancy(widx,pi)/getWardCapacity(widx));
}

double EntireSystem::rejectionProbability(int widx, vector<double> &pi){
    
    vector<double> margDist = marginalDist(widx,pi);
    
    return(margDist[margDist.size()-1]);
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