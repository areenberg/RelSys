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
 * File:   HeuristicQueue.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on September 14, 2020, 7:41 PM
 */

#include "HeuristicQueue.h"
#include "HyperQueue.h"
#include "Combinatorial.h"
#include "StatusBar.h"

#include <vector>
#include <algorithm>  
#include <iostream>
#include <limits.h>

using namespace std;

HeuristicQueue::HeuristicQueue(int main_widx,vector<vector<int>> binMap, int c, vector<int> upperLim,
        vector<int> lowerLim, double aRate, double sRate, int nhq, HyperQueue* hbQueues, QueueData * wards):
wards_pointer(wards),
hbQueues_pointer(hbQueues),
Nh(nhq),
cap(c),
arrivalRate(aRate),
serviceRate(sRate),
upperLim(upperLim),
lowerLim(lowerLim),
binMap(binMap),
main_widx(main_widx)        
{
    checkInput();
    calculateDischargeRates();
    calculateSize();
    initializeState();
    initializeJumbVectors();
}

HeuristicQueue::HeuristicQueue(const HeuristicQueue& orig) {
}

HeuristicQueue::~HeuristicQueue() {
}

void HeuristicQueue::adjustLimits(){
    //adjust limits according to the
    //bin map
    
    vector<int> ul(nBins,0);
    
    //adjust the upper limits
    for (int bidx=0; bidx<nBins; bidx++){
        for (int pidx=0; pidx<binMap.size(); pidx++){
            if (binMap[pidx][bidx]==1){
                ul[bidx] += upperLim[pidx];
            }
        }
        if (ul[bidx]>cap){
            ul[bidx]=cap;
        }
    }
    
    upperLim.resize(nBins,0);
    for (int bidx=0; bidx<nBins; bidx++){
        upperLim[bidx] = ul[bidx];
    }
    
    //adjust the lower limits
    int mn;
    lowerLim.resize(nBins,0);
    for (int bidx=0; bidx<nBins; bidx++){
        mn = numeric_limits<int>::max();
        for (int pidx=0; pidx<binMap.size(); pidx++){
            if (lowerLim[pidx]<mn){
                mn = lowerLim[pidx];
            }
        }
        lowerLim[bidx] = mn;
    }
}

void HeuristicQueue::newbinDischargeRates(vector<double> newBinDisRates){
    binDischargeRates.resize(nBins,0);
    for (int bidx=0; bidx<nBins; bidx++){
        binDischargeRates[bidx]=newBinDisRates[bidx];
    }
}

void HeuristicQueue::calculateDischargeRates(){
    //calculate the aggregated discharge rates
    
    binDischargeRates.resize(nBins,0);
    double smarr;
    for (int bidx=0; bidx<nBins; bidx++){
        smarr=0;
        for (int pidx=0; pidx<binMap.size(); pidx++){
            if (binMap[pidx][bidx]==1){
                smarr += getWardArrivalRate(pidx);  
            }
        }
        for (int pidx=0; pidx<binMap.size(); pidx++){
            if (binMap[pidx][bidx]==1){
                binDischargeRates[bidx] += (getWardArrivalRate(pidx)/smarr)*getWardServiceRate(pidx);
            }            
        }
    }
    
}

void HeuristicQueue::checkInput(){

    if (upperLim.size()!=lowerLim.size()){
        cout << "Error: Limits must be of equal length." << endl;
    }
    
    nBins = binMap[0].size();
    
    if (nBins>(Nh+1)){
        cout << "Error: There are more state bins than queues in the system." << endl;
    }
    
    adjustLimits();
    
    if (upperLim.size()!=nBins || lowerLim.size()!=nBins){
        cout << "Error: Limits must have same length as number of state bins." << endl;
    }
    int k = 0;
    for (int i=0; i<upperLim.size(); i++){
        k += lowerLim[i];
        if (upperLim[i]>cap || lowerLim[i]>cap){
            cout << "Error: Limits must be equal to or smaller than the capacity." << endl;
            
        }else if (k>cap){
            cout << "Error: Sum of lower limits must be equal to or smaller than the capacity." << endl;
        }
    }
    
}

void HeuristicQueue::marginalDist(vector<double> &pi){
    //derives the marginal distribution from the
    //overall state probability distribution.
    
    margDist.resize(margDist.size(),0);
    
    int l = 0,u = 0;
    for (int i=0; i<lowerLim.size(); i++){
        u += upperLim[i]; l += lowerLim[i];
    }
    
    int occupancy;
    initializeState();
    for (int i=0; i<Ns; i++){
        occupancy = 0;
        for (int j=0; j<lowerLim.size(); j++){
            occupancy += state[j];
        }
        margDist[occupancy-l] += pi[i];
        nextCurrentState();
    }
    
}

double HeuristicQueue::expectedOccupancy(){
    
    int l = 0,u = 0;
    for (int i=0; i<lowerLim.size(); i++){
        u += upperLim[i]; l += lowerLim[i];
    }
    
    double y = 0; int k = 0;
    for (int i=l; i<=min(u,cap); i++){
        y += i*margDist[k];
        k++;
    }
    
    return(y);
}

double HeuristicQueue::rejectionProbability(){
    return(margDist[margDist.size()-1]);
}

double HeuristicQueue::expectedOccupancyFraction(){
    return(expectedOccupancy()/cap);
}

void HeuristicQueue::buildChain(){
    //creates and stores the entire transition rate matrix
    
    qColumnIndices.clear(); qColumnIndices.resize(Ns);
    qValues.clear(); qValues.resize(Ns);
    
    for (int i=0; i<Ns; i++){
        allOutgoing();
        qColumnIndices[i].resize(toIdxSize,0);
        qValues[i].resize(toIdxSize,0);
        for (int j=0; j<toIdxSize; j++){
            qColumnIndices[i][j] = jumpToIdx[j];
            qValues[i][j] = jumpToRate[j];
        }
        nextCurrentState();
    }
    
}

void HeuristicQueue::buildTransposedChain(){
    //creates and stores the entire <transposed> transition rate matrix
    
    qColumnIndices.clear(); qColumnIndices.resize(Ns);
    qValues.clear(); qValues.resize(Ns);
    
    StatusBar sbar(Ns,30);
    for (int i=0; i<Ns; i++){
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

void HeuristicQueue::calculateSize(){
    //calculate and store the size of the state space
    
    Ns = cmb.capWithLimits(cap,upperLim,lowerLim); //size of the main queue itself;
    
    for (int i=0; i<Nh; i++){
        Ns *= getHyperSize(i);
    }
    
    //calculate size of marginal distribution
    int l=0,u=0;
    for (int i=0; i<lowerLim.size(); i++){
        u += upperLim[i]; l += lowerLim[i];
    }
    int dsize=min(u,cap)-l+1;
    margDist.resize(dsize,0);
    
    
    //hyper queue indices
    hyperWidx_vector.resize(Nh,0);
    int added = 0;
    for (int i=0; i<binMap.size(); i++){
        if (i!=main_widx){
            hyperWidx_vector[added] = i;
            added++;
        }
    }
    
    //vectors used in generating the transition rates
    fwvec.resize(nBins,0);
    hblocked.resize(Nh,false);
    
    //cout << "entire state space has size: " << Ns << " states" << endl;
    
}

void HeuristicQueue::initializeState(){
    //initialize the state
    
    sidx = 0;
    state.clear(); state.resize((lowerLim.size()+Nh),0);
    K_use = 0;
    for (int i=0; i<lowerLim.size(); i++){
        state[i] = lowerLim[i]; //elements accounting for the main queue itself
        K_use += lowerLim[i];
    }
    for (int i=lowerLim.size(); i<(lowerLim.size()+Nh); i++){
        state[i] = 0; //elements accounting for the hyper queues. starts with the blocking states
    }
    
}

void HeuristicQueue::initializeJumbVectors(){
    //initialize the jump vectors
    
    toIdxSize = 0;
    fromIdxSize = 0;
    
    int maximumSize = 1;
    int hq=Nh-1;
    for (int i=(state.size()-1); i>=lowerLim.size(); i--){
        maximumSize += max(hyperOpenStates(hq),hyperBlockedStates(hq));
        hq--;
    }
    //jumps associated with queue 
    maximumSize += lowerLim.size()*2;
    jumpToIdx.clear(); jumpToIdx.resize(maximumSize,-1);
    jumpFromIdx.clear(); jumpFromIdx.resize(maximumSize,-1);
    jumpToRate.clear(); jumpToRate.resize(maximumSize,-1);
    jumpFromRate.clear(); jumpFromRate.resize(maximumSize,-1);
    
    maximumNonZero();
    
}

void HeuristicQueue::maximumNonZero(){
    //derives an upper bound for the number of non-zero elements in the
    //transition rate matrix. This is used to assess if parameters should
    //be stored or calculated on demand.
    //This can also be used to estimate the sparsity of the matrix.
    
    maxNz = Ns * jumpToIdx.size();
    
}

void HeuristicQueue::nextCurrentState(){
    //advance to the next current state
    
    //advance state index
    if (sidx<(Ns-1)){ 
        sidx++;
    }else{
        sidx = 0;
    }
    
    //advance state form
    int i = state.size()-1;
    while (i>=0){
        
        if (i>=lowerLim.size()){ //hyper queues
            
            if (state[i] < (getHyperSize(i-lowerLim.size())-1) ){
                state[i]++;
                i = -1; //exit
            }else{
                state[i] = 0;
                i--;
            }
            
        }else{
            
            if (state[i]<upperLim[i] && K_use<cap){
                state[i]++; K_use++;
                i = -1; //exit
            }else if(state[i]>lowerLim[i]){
                K_use -= state[i]-lowerLim[i]; state[i] = lowerLim[i];
                i--;
            }else{
                i--;
            }
        }   
    }
    
}


void HeuristicQueue::previousCurrentState(){
    //go back to the previous current state
    
    //advance state index
    if (sidx>0){ 
        sidx--;
    }else{
        sidx = Ns-1;
    }
    
    //advance state form
    int i = state.size()-1;
    while (i>=0){
        
        if (i>=lowerLim.size()){ //hyper queues
            
            if (state[i]>0){
                state[i]--;
                i = -1; //exit
            }else{
                state[i] = getHyperSize(i-lowerLim.size())-1;
                i--;
            }
            
        }else{
            
            if (state[i]>lowerLim[i] && K_use>0  <upperLim[i] && K_use<cap){
                state[i]--; K_use--;
                i = -1; //exit
            }else if(state[i]<upperLim[i]){
                K_use += upperLim[i]-state[i]; state[i] = upperLim[i];
                i--;
            }else{
                i--;
            }
        }   
    }
    
}


void HeuristicQueue::allOutgoing(){
    //all hyper queues can jump to an open/blocked state.
    //queue nodes can discharge or admit a customer when:
    //1. They do not violate the capacity or user-specified limits.
    //2. The nodes linked to the hyper queues can only admit when the
    //associated hyper queue is in a blocked state.
    
    int hq, jump, prod;
    double diag = 0;
    bool hasMain,jmpallow;
    toIdxSize = 0;
    
    //hyper queue jumps
    hq=Nh-1; prod = 1;
    for (int i=(state.size()-1); i>=lowerLim.size(); i--){
        
        if(state[i]<hyperBlockedStates(hq)){ //state is blocked
            hblocked[hq]=true;
            for (int j=0; j<hyperOpenStates(hq); j++){ //loop over all open states the process can jump to
                jump = (j+hyperBlockedStates(hq)) - state[i];
                jump *= prod;
                jumpToIdx[toIdxSize] = sidx+jump;
                
                jumpToRate[toIdxSize] = getHyperBlockedRate(hq,state[i])*getHyperOpenDist(hq,j);
                diag -= jumpToRate[toIdxSize];
                
                toIdxSize++; 
            }
            
        }else{ //state is open
            hblocked[hq]=false;
            for (int j=0; j<hyperBlockedStates(hq); j++){ //loop over all blocked states the process can jump to
                jump = j + hyperOpenStates(hq) - (getHyperSize(hq)-1-state[i]); 
                jump *= prod;
                jumpToIdx[toIdxSize] = sidx-jump;
                
                jumpToRate[toIdxSize] = getHyperOpenRate(hq, (state[i]-hyperBlockedStates(hq)) )*getHyperBlockedDist(hq,j);
                diag -= jumpToRate[toIdxSize];
                
                toIdxSize++;
            }
        }
        
        prod *= getHyperSize(hq);
        hq--;
    }
    
    //queue jumps
    //admission
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
        //ihq = i+Nh;
        //hq = ihq - lowerLim.size();
        
        //check if bin is associated with main queue or has blocked hyper queue
        jmpallow=false;
        hasMain=false;
        hq=0;
        if (binMap[main_widx][bidx]==1){
            hasMain=true;
            jmpallow=true;
        }
        while (jmpallow==false && hq<hyperWidx_vector.size()){
            if (binMap[hyperWidx_vector[hq]][bidx]==1 && hblocked[hq]){
                jmpallow=true;
            }
            hq++;
        }
        
        //if (bidx==0 || state[ihq]<hyperBlockedStates(hq)  ){
        if (jmpallow && state[bidx]<upperLim[bidx] && K_use<cap){
            //if (state[bidx]<upperLim[bidx] && K_use<cap){
                //admit one
                if (bidx<(lowerLim.size()-1)){
                    jump = forwardOne(K_use,fwvec,(state[bidx]+1),bidx);
                }else{
                    jump = 1;
                }
                jump *= prod;
                jumpToIdx[toIdxSize] = sidx+jump;
                
//                if (i==0){
//                    jumpToRate[toIdxSize] = arrivalRate;
//                }else{
//                    jumpToRate[toIdxSize] = getHyperArrivalRate(hq);
//                }
                
                jumpToRate[toIdxSize] = getAdmissionRate(bidx,hblocked,hasMain);
                
                diag -= jumpToRate[toIdxSize];
                
                toIdxSize++;
           // }
        }
        
    }
    //discharge
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
        if (state[bidx]>lowerLim[bidx]){
            //discharge one
            if (bidx<(lowerLim.size()-1)){
                jump = backwardOne(K_use,fwvec,(state[bidx]-1),bidx);
            }else{
                jump = 1;
            }
            jump *= prod;
            jumpToIdx[toIdxSize] = sidx-jump;
            
//            if (bidx==0){
//                jumpToRate[toIdxSize] = serviceRate*state[bidx];
//            }else{
//                hq = bidx-1;
//                jumpToRate[toIdxSize] = getHyperServiceRate(hq)*state[bidx]; 
//            }
            
            jumpToRate[toIdxSize] = binDischargeRates[bidx]*state[bidx];
            
            diag -= jumpToRate[toIdxSize];
            
            toIdxSize++;
        }
    }
    //insert diagonal
    jumpToIdx[toIdxSize] = sidx;
    jumpToRate[toIdxSize] = diag;
    toIdxSize++;
    
    
}

void HeuristicQueue::allIngoing(){
    //all hyper queues can come from an open/blocked state.
    //queue nodes can come from a discharge or admission of a customer when:
    //1. The old state did not violate the capacity or user-specified limits.
    //2. The nodes linked to the hyper queues can only admit when the
    //associated hyper queue is in a blocked state.
    
    int hq, jump, prod; //delta, l, jj, k;
    bool hasMain,jmpallow;
    fromIdxSize = 0;
    
    //hyper queue jumps
    hq=Nh-1; prod = 1;
    hblocked.resize(Nh,false);
    for (int i=(state.size()-1); i>=lowerLim.size(); i--){
        
        if(state[i]<hyperBlockedStates(hq)){ //state is blocked
            hblocked[hq]=true;
            for (int j=0; j<hyperOpenStates(hq); j++){ //loop over all open states the process <came from>
                jump = (j+hyperBlockedStates(hq)) - state[i];
                jump *= prod;
                jumpFromIdx[fromIdxSize] = sidx+jump;
                //cout << "from open: " << jumpFromIdx[fromIdxSize] << endl;
                
                jumpFromRate[fromIdxSize] = getHyperOpenRate(hq,j)*getHyperBlockedDist(hq,state[i]);
                
                fromIdxSize++; 
            }
            
        }else{ //state is open
            hblocked[hq]=false;
            for (int j=0; j<hyperBlockedStates(hq); j++){ //loop over all blocked states the process <came from>
                jump = j + hyperOpenStates(hq) - (getHyperSize(hq)-1-state[i]); 
                jump *= prod; 
                jumpFromIdx[fromIdxSize] = sidx-jump;
                //cout << "from blocked: " << jumpFromIdx[fromIdxSize] << endl;
                
                jumpFromRate[fromIdxSize] = getHyperBlockedRate(hq,j)*getHyperOpenDist(hq, (state[i]-hyperBlockedStates(hq)) );
                
                fromIdxSize++;
            }
        }
        
        prod *= getHyperSize(hq);
        hq--;
    }
    
    //queue jumps    
    //came from a discharge
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
         
        if (state[bidx]<upperLim[bidx] && K_use<cap){
            //find jump size up to the previous state
            if (bidx<(lowerLim.size()-1)){
                jump = forwardOne(K_use,fwvec,(state[bidx]+1),bidx);                
            }else{
                jump = 1;
            }
            
            jump *= prod;
            jumpFromIdx[fromIdxSize] = sidx+jump;
//            cout << "from discharge: " << jumpFromIdx[fromIdxSize] << endl;
//            cout << "sidx=" << sidx << ", jump=" << jump << endl;
            
//            if (bidx==0){
//                jumpFromRate[fromIdxSize] = serviceRate*(state[bidx]+1);
//            }else{
//                hq = bidx-1;
//                jumpFromRate[fromIdxSize] = getHyperServiceRate(hq)*(state[bidx]+1);
//            }
            
            jumpFromRate[fromIdxSize] = binDischargeRates[bidx]*(state[bidx]+1);
                
            fromIdxSize++;
        }
        
    }
    
    
    
    //came from an admission
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
//        ihq = bidx+Nh;
//        hq = ihq - lowerLim.size();
        
        //check if bin is associated with main queue or has blocked hyper queue
        jmpallow=false;
        hasMain=false;
        hq=0;
        if (binMap[main_widx][bidx]==1){
            hasMain=true;
            jmpallow=true;
        }
        while (jmpallow==false && hq<hyperWidx_vector.size()){
            if (binMap[hyperWidx_vector[hq]][bidx]==1 && hblocked[hq]){
                jmpallow=true;
            }
            hq++;
        }
        
        if (jmpallow && state[bidx]>lowerLim[bidx]){
//        if (state[bidx]>lowerLim[bidx] && (bidx==0 || state[ihq]<hyperBlockedStates(hq)) ){
            //find jump size down to the previous state
            if (bidx<(lowerLim.size()-1)){
                jump = backwardOne(K_use,fwvec,(state[bidx]-1),bidx);
            }else{
                jump = 1;
            }
            jump *= prod;
            jumpFromIdx[fromIdxSize] = sidx-jump;
            //cout << "from admission: " << jumpFromIdx[fromIdxSize] << endl;
            
//            if (bidx==0){
//                jumpFromRate[fromIdxSize] = arrivalRate;
//            }else{
//                hq = bidx-1;
//                jumpFromRate[fromIdxSize] = getHyperArrivalRate(hq); 
//            }
            
            jumpFromRate[fromIdxSize] = getAdmissionRate(bidx,hblocked,hasMain);
            
            fromIdxSize++;
        }
    }
    
    //calculate and insert diagonal
    jumpFromIdx[fromIdxSize] = sidx;
    jumpFromRate[fromIdxSize] = calculateDiagonal();
    fromIdxSize++;
    
}

double HeuristicQueue::getAdmissionRate(int &bidx, vector<bool> &hblocked, bool &hasMain){
    //returns the rate with which admissions
    //occur in the specified bin index (bidx).
    
    double rate;
    if (hasMain){
        rate=arrivalRate;
    }else{
        rate=0;
    }
    
    for (int hq=0; hq<Nh; hq++){
        if (binMap[hyperWidx_vector[hq]][bidx]==1 && hblocked[hq]){
            rate+=getHyperArrivalRate(hq);
        }
    }
    
    return(rate);
}

double HeuristicQueue::calculateDiagonal(){
    //calculates the value of the diagonal related
    //to the current state.
    
    double diag=0;
    
//    int hq, ihq;
    
//    //hyper queue jumps
//    hq=Nh-1;
//    for (int i=(state.size()-1); i>=lowerLim.size(); i--){
//        
//        if(state[i]<hyperBlockedStates(hq)){ //state is blocked
//            for (int j=0; j<hyperOpenStates(hq); j++){ //loop over all open states the process can jump to
//                diag -= getHyperBlockedRate(hq,state[i])*getHyperOpenDist(hq,j); 
//            }
//            
//        }else{ //state is open
//            for (int j=0; j<hyperBlockedStates(hq); j++){ //loop over all blocked states the process can jump to
//                diag -= getHyperOpenRate(hq, (state[i]-hyperBlockedStates(hq)) )*getHyperBlockedDist(hq,j);
//            }
//        }
//        
//        hq--;
//    }
//    
//    //queue jumps
//    //admission
//    for (int i=(lowerLim.size()-1); i>=0; i--){
//        ihq = i+Nh;
//        hq = ihq - lowerLim.size(); 
//        if (i==0 || state[ihq]<hyperBlockedStates(hq)  ){
//            if (state[i]<upperLim[i] && K_use<cap){
//                
//                if (i==0){
//                    diag -= arrivalRate;
//                }else{
//                    diag -= getHyperArrivalRate(hq);
//                }
//                
//            }
//        }
//        
//    }
//    //discharge
//    for (int i=(lowerLim.size()-1); i>=0; i--){
//        if (state[i]>lowerLim[i]){
//            
//            if (i==0){
//                diag -= serviceRate*state[i];
//            }else{
//                hq = i-1;
//                diag -= getHyperServiceRate(hq)*state[i]; 
//            }
//            
//        }
//    }
    
    
    int hq;
    bool hasMain,jmpallow;
    
    //hyper queue jumps
    hq=Nh-1;
    for (int i=(state.size()-1); i>=lowerLim.size(); i--){
        
        if(state[i]<hyperBlockedStates(hq)){ //state is blocked
            hblocked[hq]=true;
            for (int j=0; j<hyperOpenStates(hq); j++){ //loop over all open states the process can jump to
                diag -= getHyperBlockedRate(hq,state[i])*getHyperOpenDist(hq,j); 
            }
            
        }else{ //state is open
            hblocked[hq]=false;
            for (int j=0; j<hyperBlockedStates(hq); j++){ //loop over all blocked states the process can jump to
                diag -= getHyperOpenRate(hq, (state[i]-hyperBlockedStates(hq)) )*getHyperBlockedDist(hq,j);
            }
        }
        
        hq--;
    }
    
    //queue jumps
    //admission
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
        
        //check if bin is associated with main queue or has blocked hyper queue
        jmpallow=false;
        hasMain=false;
        hq=0;
        if (binMap[main_widx][bidx]==1){
            hasMain=true;
            jmpallow=true;
        }
        while (jmpallow==false && hq<hyperWidx_vector.size()){
            if (binMap[hyperWidx_vector[hq]][bidx]==1 && hblocked[hq]){
                jmpallow=true;
            }
            hq++;
        }
        
        if (jmpallow && state[bidx]<upperLim[bidx] && K_use<cap){
                diag -= getAdmissionRate(bidx,hblocked,hasMain);
        }
        
    }
    //discharge
    for (int bidx=(lowerLim.size()-1); bidx>=0; bidx--){
        if (state[bidx]>lowerLim[bidx]){
            //discharge one
            diag -= binDischargeRates[bidx]*state[bidx];
        }
    }
    
    return(diag);
}



int HeuristicQueue::forwardOne(int Ku, vector<int> &j, int targetval, int targetidx){
    
    for (int i=0; i<lowerLim.size(); i++){
        j[i] = state[i];
    }
    
    int i, ii, m = 0;
    
    bool run = true;
    while(run){
        
        //advance state form by one step
        m++;
        
        i = lowerLim.size()-1;
        while (i>=0){
            if (j[i]<upperLim[i] && Ku<cap){
                j[i]++; Ku++;
                i = -1; //exit
            }else if(j[i]>lowerLim[i]){
                Ku -= j[i]-lowerLim[i]; j[i] = lowerLim[i];
                i--;
            }else{
                i--;
            }
        }
        
        if ( j[targetidx]==targetval ){
            ii=0; run = false;
            do{   
                if (j[ii]!=state[ii] && ii!=targetidx){
                    run = true;
                }
                ii++;
            }while(run==false && ii<lowerLim.size());
            
        }
        
    }
    
    return (m);
}


int HeuristicQueue::backwardOne(int Ku, vector<int> &j, int targetval, int targetidx){
    
    for (int i=0; i<lowerLim.size(); i++){
        j[i] = state[i];
    }
    j[targetidx] = targetval;
    Ku--;
    
    int i, ii, m = 0;
    
    bool run = true;
    while(run){
        
        //advance state form by one step
        m++;
        
        i = lowerLim.size()-1;
        while (i>=0){
            if (j[i]<upperLim[i] && Ku<cap){
                j[i]++; Ku++;
                i = -1; //exit
            }else if(j[i]>lowerLim[i]){
                Ku -= j[i]-lowerLim[i]; j[i] = lowerLim[i];
                i--;
            }else{
                i--;
            }
        }
        
        if ( j[targetidx]==state[targetidx] ){
            ii=0; run = false;
            do{   
                if (j[ii]!=state[ii]){
                    run = true;
                }
                ii++;
            }while(run==false && ii<lowerLim.size());
            
        }
        
    }
    
    return (m);
}



int HeuristicQueue::hyperOpenStates(int hq){
    
    return((hbQueues_pointer + hq)->openRates.size());
}

int HeuristicQueue::hyperBlockedStates(int hq){
    
    return((hbQueues_pointer + hq)->blockedRates.size());
}

double HeuristicQueue::getHyperOpenRate(int hq, int idx){
    
    return((hbQueues_pointer + hq)->openRates[idx]);
}

double HeuristicQueue::getHyperOpenDist(int hq, int idx){
    
    return((hbQueues_pointer + hq)->openDist[idx]);
}

double HeuristicQueue::getHyperBlockedRate(int hq, int idx){
    
    return((hbQueues_pointer + hq)->blockedRates[idx]);
}

double HeuristicQueue::getHyperBlockedDist(int hq, int idx){
    
    return((hbQueues_pointer + hq)->blockedDist[idx]);
}

int HeuristicQueue::getHyperSize(int hq){
    
    return((hbQueues_pointer + hq)->numberOfStates);
}

double HeuristicQueue::getHyperArrivalRate(int hq){
    
    return((hbQueues_pointer + hq)->arrivalRate);
}

double HeuristicQueue::getHyperServiceRate(int hq){
    
    return((hbQueues_pointer + hq)->serviceRate);
}

double HeuristicQueue::getWardArrivalRate(int ward){
    
    return((wards_pointer + ward)->arrivalRate);
}

double HeuristicQueue::getWardServiceRate(int ward){
    
    return((wards_pointer + ward)->serviceRate);
}