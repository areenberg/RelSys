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
 * File:   RelocEvaluation.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on November 24, 2020, 9:31 PM
 */

#include "RelocEvaluation.h"
#include "HeuristicQueue.h"
#include "HyperQueue.h"
#include "LinSolver.h"
#include "QueueData.h"
#include "RelocSimulation.h"
#include "PhaseFitter.h"

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

RelocEvaluation::RelocEvaluation(int nW, QueueData * wards):
wards_pointer(wards),
nWards(nW),
simReady(false),
validateReady(false),
simMargDist(false),
changeDisRates(false)        
{
    initializeSystem();
}

RelocEvaluation::RelocEvaluation(const RelocEvaluation& orig) {
}

RelocEvaluation::~RelocEvaluation() {
}

void RelocEvaluation::initializeSystem(){

    //set states in hyper queues
    setDefaultHyperQueueStates();
    
    //simulation variables
    sim_pointer = new RelocSimulation(nWards,wards_pointer);
    
    //set default binMap
    setDefaultBinMap();
    
}

void RelocEvaluation::setDefaultHyperQueueStates(){
    hyperOpenStates.resize((nWards-1),2); //2 open states
    hyperBlockedStates.resize((nWards-1),1); //1 blocked state
}

void RelocEvaluation::setOpenHyperStates(int n){
    //modifies all open hyper queues
    for (int i=0; i<hyperOpenStates.size(); i++){
        hyperOpenStates[i] = n;
    }
}

void RelocEvaluation::setBlockedHyperStates(int n){
    //modifies all blocked hyper queues
    for (int i=0; i<hyperBlockedStates.size(); i++){
        hyperBlockedStates[i] = n;
    }    
}

void RelocEvaluation::setDefaultBinMap(){
    //automatically merge patient types
    //with equal service rates
    
    vector<bool> unusedColumn(nWards,true);
    int idx,j;
    bool check;
    
    binMap.resize(nWards);
    for (int i=0; i<nWards; i++){
        idx=i;
        check=true;
        binMap[i].resize(nWards,0);
        j=0;
        while (check && j<nWards){
            if (check && i!=j && getWardServiceRate(i)==getWardServiceRate(j) && j<i){
                idx=j;
                check=false;
            }
            j++;
        }
        
        binMap[i][idx]=1;
        unusedColumn[idx]=false;
    }
    
    int inUse=0;
    for (int i=0; i<nWards; i++){
        if (unusedColumn[i]==false){
            inUse++;
        }
    }
    
    int k;
    for (int i=0; i<nWards; i++){
        binMap[i].resize(inUse,0);
        k=0;
        for (int j=0; j<nWards; j++){
            if (unusedColumn[j]==false){
                binMap[i][k]=binMap[i][j];
                k++;
            }
            
        }
    }
    //cout << "State reduced with " << (nWards-binMap[0].size()) << " bins." << endl;

}


void RelocEvaluation::setBinMap(vector<vector<int>> bMap){
    //overwrite the default binMap
    
    binMap.resize(bMap.size());
    for (int widx=0; widx<bMap.size(); widx++){
        binMap[widx].resize(bMap[widx].size(),0);
        for (int bidx=0; bidx<bMap[widx].size(); bidx++){
            binMap[widx][bidx]=bMap[widx][bidx];
        }
    }
    
}

void RelocEvaluation::setBinDischargeRates(vector<double> disRates){
    changeDisRates=true;
    newDisRates.resize(disRates.size(),0);
    for (int i=0; i<disRates.size(); i++){
        newDisRates[i]=disRates[i];
    }
}


void RelocEvaluation::runSimulation(int sd, int burnIn,
        int minTime, int minSamples){
    
    seed = sd;
    
    sim_pointer->setSeed(seed);
    
    //auto start = high_resolution_clock::now();
    
    vector<int> maxWardSamples(nWards,100000);
    sim_pointer->simulate(burnIn,minTime,maxWardSamples,minSamples);
    
    //auto stop = high_resolution_clock::now(); //stop time 
    //auto duration = duration_cast<milliseconds>(stop - start); 
    //cout << "Runtime of simulation: " << duration.count() << " milliseconds\n" << endl;
    
    simReady = true;
}

bool RelocEvaluation::isSingleWard(int widx){
    
    if (nWards==1||noRelocations(widx)){
        return(true);
    }else{
        return(false);
    }
    
}


bool RelocEvaluation::noRelocations(int widx){
    
    bool noReloc=true;
    for (int j=0; j<nWards; j++){
        if (widx!=j && getWardRelocationProbabilities(j)[widx]>0.0 && getWardArrivalRate(j)>1e-16){
            noReloc=false;
        }
    }
    return(noReloc);
}


void RelocEvaluation::evalSingleWard(int widx){
    marginalDist.resize((getWardCapacity(widx)+1),0);
    marginalDistPref.resize((getWardCapacity(widx)+1),0);
    if (getWardArrivalRate(widx)<=1e-16){
        marginalDist[0] = 1.0;
        marginalDistPref[0] = 1.0;
        blockingProbability = 0.0;
        blockingProbabilityPref = 0.0;
        
        expectedOccupancy = 0.0;
        expOccFraction = 0.0;
    }else{
        //use the erlang loss model
        double lambda = getWardArrivalRate(widx);
        double mu = getWardServiceRate(widx);
        int servers = getWardCapacity(widx);
        for (int k=0; k<marginalDist.size(); k++){
            marginalDist[k] = erlangLoss(k,lambda,mu,servers);
            marginalDistPref[k] = marginalDist[k];
        }
        
        expectedOccupancy = 0.0;
        for (int k=0; k<marginalDistPref.size(); k++){
            expectedOccupancy += marginalDistPref[k]*k;
        }
        expOccFraction = expectedOccupancy/(double)servers;
        blockingProbability = marginalDist[marginalDist.size()-1];
        blockingProbabilityPref = marginalDistPref[marginalDistPref.size()-1];
    }
}

void RelocEvaluation::validateModel(int main_widx){
    
    mfocus = main_widx;
    
    if (isSingleWard(mfocus)){
        
        stateSpaceSize = getWardCapacity(mfocus);
        validateReady = true; //tag the model as validated
        
    }else if (simReady){
    
        //create and add surrogate queues to the system
        double arr;
        nhq = nWards-1; //number of hyper queues
        
        //upper and lower limits in each queue. must include all queues (main and hyper queues)
        upperLimits.resize(nWards,0);
        lowerLimits.resize(nWards,0);
        //automatically adjust truncation
        //cout << "Setting limits... ";
        setUpperLimits(upperLimits,mfocus);
//        setLowerLimits(lowerLimits,mfocus);
        //cout << "done." << endl;
        
        //hyper queue indices
        hyperWidx_vector.resize(nhq,0);
        int added = 0;
        for (int i=0; i<nWards; i++){
            if (i!=mfocus){
                hyperWidx_vector[added] = i;
                added++;
            }
        }
        
        //create objects for each surrogate (hyper) queue
        hq_array = new HyperQueue[nhq];
        for (int i=0; i<nhq; i++){
            arr =  getWardArrivalRate(hyperWidx_vector[i])*getWardRelocationProbabilities(hyperWidx_vector[i])[mfocus];
            
            hq_array[i] = HyperQueue(hyperWidx_vector[i],hyperBlockedStates[i],hyperOpenStates[i],//statesBlocked,statesOpen,
                arr,getWardServiceRate(hyperWidx_vector[i]),sim_pointer);
        }
        
        //create the main (heuristic) queue object
        hqueue = new HeuristicQueue(mfocus,binMap,getWardCapacity(mfocus),upperLimits,lowerLimits,
                getWardArrivalRate(mfocus),getWardServiceRate(mfocus),nhq,hq_array,wards_pointer);
        if (changeDisRates){
            //change the bin discharge rates
            hqueue->newbinDischargeRates(newDisRates);
        }
        
        cout << "Number of states = " << hqueue->Ns << endl;
        stateSpaceSize = hqueue->Ns;
        
        validateReady = true; //tag the model as validated
    }else{
        cout << "Execution failed because simulation was not available." << endl;
        cout << "In order to solve this issue, run the runSimulation() method prior to runHeuristic()." << endl;
    }    
        
}


void RelocEvaluation::evaluateModel(){
    if (isSingleWard(mfocus)){
        
        cout << "Detected single queue - evaluate using Erlang B." << endl;
        evalSingleWard(mfocus);
        
    }else if (validateReady){
        
        //fit PH parameters for each surrogate (hyper) queue
        cout << "Fitting PH parameters... " << endl;
        for (int i=0; i<nhq; i++){
            cout << "Fit " << (i+1) << endl;
            hq_array[i].fitAll(seed);
        }
        cout << "Done." << endl;
        
        //cout << "Evaluating Queue " << (mfocus+1) << "..." << endl;
        
        //solve for steady-state distribution
        LinSolver solver;
        initializeStateDistribution(hqueue);    
        double solverTolerance = 1e-9;
        solver.sor(pi,hqueue,1.0,solverTolerance);
            
        vmemory = solver.vmemory; //estimate of maximum memory usage
            
        //store the marginal distribution and some other metrics
        hqueue->marginalDist(pi); //calculate both marginal dist. types
        marginalDistPref = hqueue->margDistPref;
        marginalDist = hqueue->margDist;
            
        expectedOccupancy = hqueue->expectedOccupancy();
        expOccFraction = expectedOccupancy/(double)getWardCapacity(mfocus);
        blockingProbability = marginalDist[marginalDist.size()-1];
        blockingProbabilityPref = marginalDistPref[marginalDistPref.size()-1];
        
        delete hqueue; //free memory related to model
        
    }else{
        cout << "Execution failed because model has not been validated." << endl;
        cout << "In order to solve this issue, run the validateModel() method prior to evaluateModel()." << endl;
    }
    
}


void RelocEvaluation::runHeuristic(int main_widx){
    
        auto start = high_resolution_clock::now(); //start time 
        
        validateModel(main_widx);
        
        //cout << "Evaluating Queue " << (main_widx+1) << "..." << endl;
        evaluateModel();
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        cout << "Runtime of CTMC (excl. simulation): " << duration.count() << " milliseconds\n" << endl;
    
}

void RelocEvaluation::setUpperLimits(vector<int> &upperLimits, int &main_widx){
    //finds the upper truncation limit for each patient
    //type in the main ward.
    //finds the smallest truncation that has at least tailden probability
    //mass in the upper tail of the distribution.
    
    //tail cut-off probability mass
    double tailden = 1e-6;
    
    vector<vector<int>> fd = sim_pointer->freqDist[main_widx];
    
    double mn,k,p;
    
    //cout << "Upper truncation cap. limits:" << endl;
    for (int pidx=0; pidx<nWards; pidx++){
        
        //truncation using Chebyshev's inequality
        mn = sampleMean(fd[pidx]);
        k = ceil(mn);
        if ((k+1)<getWardCapacity(main_widx)){
            do{
                k++;
                p = chebyshevBound(k,fd[pidx]);
            }while(k<getWardCapacity(main_widx) && p>tailden);
            upperLimits[pidx] = k;
        }else{
            upperLimits[pidx] = getWardCapacity(main_widx);
        }
            
        //cout << upperLimits[pidx] << " " << flush;
    }
    //cout << endl;
    
}

//void RelocEvaluation::setLowerLimits(vector<int> &lowerLimits, int &main_widx){
//    //lower truncation limits are not adjusted (for now).
//    
//    //cout << "Lower truncation cap. limits:" << endl;
//    for (int pidx=0; pidx<nWards; pidx++){
//        //cout << lowerLimits[pidx] << " " << flush;
//    }
//    //cout << endl;
//    
//}

double RelocEvaluation::sampleMean(vector<int> &freqDist){
    
    double y,sm;
    sm=0;
    for (int i=0; i<freqDist.size(); i++){
        sm += freqDist[i];
    }
    y=0;
    for (int i=0; i<freqDist.size(); i++){
        y += (freqDist[i]/sm)*i;
    }
    
    return(y);
}

double RelocEvaluation::sampleSD(vector<int> &freqDist){
    
    double n, xhat, sqdiff, diff;
    xhat = sampleMean(freqDist);
    n=0;
    for (int i=0; i<freqDist.size(); i++){
        n += freqDist[i];
    }
    sqdiff = 0;
    for (int i=0; i<freqDist.size(); i++){
        for (int j=0; j<freqDist[i]; j++){
            diff = i-xhat;
            sqdiff += pow(diff,2.0);
        }
    }
    
    return(sqrt((1.0/(n-1.0))*sqdiff));
}

double RelocEvaluation::chebyshevBound(double &x, vector<int> &freqDist){
    //sample version of Chebyshev's inequality. 
    //assumes input variable x is larger than the mean
    //of freqDist.
    
    double n,p,mn,sd,k,ginput;
    mn = sampleMean(freqDist);
    sd = sampleSD(freqDist);
    n=0;
    for (int i=0; i<freqDist.size(); i++){
        n += freqDist[i];
    }
    k = (x-mn)/sd;
    ginput = (n*pow(k,2.0))/(n-1+pow(k,2.0));
    
    p = Gfunction((n+1),ginput)/(n+1);
    p *= pow((n/(n+1)),0.5);
    
    return(p);
}

double RelocEvaluation::Gfunction(double q, double x){
    
    double a, R;
    R = floor(q/x);
    a = (q*(q-R))/(1.0+R*(q-R));
    
    if ((int)R%2==0){
        return(R); 
    }else if(x<a){
        return(R);
    }else{
        return(R-1);
    }
    
}

void RelocEvaluation::initializeStateDistribution(HeuristicQueue * hqueue, bool erlangInit){
    //initial state distribution
    pi.resize(hqueue->Ns,0);
    
    double sm=0;
    if (erlangInit){
        //initialize using the Erlang loss model
        hqueue->initializeState();
        int K;
        for (int i=0; i<hqueue->Ns; i++){
            K=0;
            for (int j=0; j<hqueue->nBins; j++){
                K+=hqueue->state[j];
            }
            pi[i]=erlangLoss(K,hqueue->arrivalRate,hqueue->serviceRate,hqueue->cap);
            sm+=pi[i];
            hqueue->nextCurrentState();
        }
    }else{
        //initialize using random numbers
        mt19937 rgen(seed);
        uniform_real_distribution<> dis(1,hqueue->Ns);
    
        for (int i=0; i<hqueue->Ns; i++){
            pi[i] = dis(rgen);
            sm += pi[i]; 
        }
    }     
    
    //normalize    
    for (int i=0; i<hqueue->Ns; i++){
        pi[i] /= sm;
    }
    
}

double RelocEvaluation::erlangLoss(int &k, double &lambda, double &mu, int &servers){
    //calculates the probability of k occupied servers
    //using the Erlang loss model
    
    double r=lambda/mu;
    
    double fr1 = pow(r,(double)k)/factorial(k);
    double fr2=0;
    for (int i=0; i<=servers; i++){
        fr2+=pow(r,(double)i)/factorial(i);
    }
    return(fr1/fr2);
}

long double RelocEvaluation::factorial(int &x){
    long double fact = 1;
    if (x>0){
        for (int i=1; i<=x; i++){
            fact *= i;
        }
    }
    return(fact);
}

int RelocEvaluation::getWardID(int ward){
    
    return((wards_pointer + ward)->wardnr);
}

double RelocEvaluation::getWardArrivalRate(int ward){
    
    return((wards_pointer + ward)->arrivalRate);
}

double RelocEvaluation::getWardServiceRate(int ward){
    
    return((wards_pointer + ward)->serviceRate);
}

vector<double> RelocEvaluation::getWardRelocationProbabilities(int ward){
    
    return((wards_pointer + ward)->relocationProbabilities);
}

//void RelocEvaluation::calculateWardStateSpaceSize(int ward, int numberOfWards){
//    
//    return((wards_pointer + ward)->calculateWardStateSpace(numberOfWards));
//}

int RelocEvaluation::getWardStateSpaceSize(int ward){
    
    return((wards_pointer + ward)->wardStateSpaceSize);
}

int RelocEvaluation::getWardCapacity(int ward){
    
    return((wards_pointer + ward)->capacity);
}