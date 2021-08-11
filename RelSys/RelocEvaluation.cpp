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
simReady(false)        
{
    initializeSystem();
}

RelocEvaluation::RelocEvaluation(const RelocEvaluation& orig) {
}

RelocEvaluation::~RelocEvaluation() {
}

void RelocEvaluation::initializeSystem(){
    
    //simulation variables
    sim_pointer = new RelocSimulation[1];
    sim_pointer[0] = RelocSimulation(nWards,wards_pointer);
    
}

void RelocEvaluation::runSimulation(int sd, int burnIn,
        int minTime, int minSamples){
    
    seed = sd;
    
    sim_pointer[0].setSeed(seed);
    
    auto start = high_resolution_clock::now();
    
    vector<int> dummyVec(1,-1);
    sim_pointer[0].simulate(burnIn,minTime,dummyVec,minSamples);
    
    auto stop = high_resolution_clock::now(); //stop time 
    auto duration = duration_cast<milliseconds>(stop - start); 
    cout << "Runtime of simulation: " << duration.count() << " milliseconds\n" << endl;
    
    simReady = true;
}

void RelocEvaluation::runHeuristic(int main_widx){
    
    if (simReady){
        auto start = high_resolution_clock::now(); //start time 
        
        //create and add surrogate queues to the system
        double arr;
        int nhq = nWards-1; //number of hyper queues
        int statesBlocked = 1;
        int statesOpen = 2;
    
        //hyper queue indices
        vector<int> hyperWidx_vector(nhq,0);
        int added = 0;
        for (int i=0; i<nWards; i++){
            if (i!=main_widx){
                hyperWidx_vector[added] = i;
                added++;
            }
        }
        
        //prepare and fit PH parameters for each surrogate (hyper) queue
        cout << "Fitting PH parameters ..." << endl;
        HyperQueue * hq_array = new HyperQueue[nhq];
        for (int i=0; i<nhq; i++){
            arr =  getWardArrivalRate(hyperWidx_vector[i])*getWardRelocationProbabilities(hyperWidx_vector[i])[main_widx];
            
            hq_array[i] = HyperQueue(hyperWidx_vector[i],statesBlocked,statesOpen,
                arr,getWardServiceRate(hyperWidx_vector[i]),sim_pointer);
            cout << "Fit " << (i+1) << endl;
            hq_array[i].fitAll(seed);
        }
        cout << "Done." << endl;
        
        //specifications of main queue
        double arrivalRate = getWardArrivalRate(main_widx);
        double serviceRate = getWardServiceRate(main_widx);
        int capacity = getWardCapacity(main_widx);
    
        //upper and lower limits in each queue. must include all queues (main and hyper queues)
        vector<int> upperLimits(nWards,0);
        vector<int> lowerLimits(nWards,0);
        //automatically adjust truncation
        cout << "Setting limits..." << endl;
        setUpperLimits(upperLimits,main_widx);
        setLowerLimits(lowerLimits,main_widx);
        
        //create the main (heuristic) queue object
        HeuristicQueue hqueue(capacity,upperLimits,lowerLimits,arrivalRate,serviceRate,nhq,hq_array);
        
        cout << "Evaluating Ward " << (main_widx+1) << "..." << endl;
        cout << "Number of states = " << hqueue.Ns << endl;
    
        //solve for steady-state distribution
        initializeStateDistribution(hqueue);
        LinSolver solver;
        double solverTolerance = 1e-9;
        solver.sor(pi,hqueue,1.0,solverTolerance);
        
        vmemory = solver.vmemory; //estimate of maximum memory usage
        
        //store the marginal distribution and some other metrics
        marginalDist = hqueue.marginalDist(pi);
        expectedOccupancy = hqueue.expectedOccupancy(pi);
        expOccFraction = hqueue.expectedOccupancyFraction(pi);
        blockingProbability = marginalDist[marginalDist.size()-1];
        
        
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        cout << "Runtime of CTMC (excl. simulation): " << duration.count() << " milliseconds\n" << endl;
        
        //print and store some general metrics
        //cout << "\nWARD " << (main_widx+1) << " STATS:" << endl;
        
        //cout << "Expected load = " << expectedOccupancy << endl;
        //cout << "Capacity utilization = " << expOccFraction*100 << "%\n" << endl;
        
        
    }else{
        cout << "Execution failed because simulation was not available." << endl;
        cout << "In order to solve this issue, run the runSimulation() method prior to runHeuristic()." << endl;
    }
    
}

void RelocEvaluation::setUpperLimits(vector<int> &upperLimits, int &main_widx){
    //finds the upper truncation limit for each patient
    //type in the main ward.
    //finds the smallest truncation that has at least tailden probability
    //mass in the upper tail of the distribution.
    
    //tail cut-off probability mass
    double tailden = 1e-3;
    
    //vector<vector<double>> dd = sim_pointer[0].denDist[main_widx];
    vector<vector<int>> fd = sim_pointer[0].freqDist[main_widx];
    
//    cout << "Freq. distribution:" << endl;
//    for (int i=0; i<fd.size(); i++){
//        for (int j=0; j<fd[i].size(); j++){
//            cout << fd[i][j] << " " << flush;
//        }
//        cout << endl;
//    }
    
    double mn,k,p;
    
    cout << "Upper truncation cap. limits" << endl;
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
            
        cout << upperLimits[pidx] << " " << flush;
    }
    cout << endl;
    
}

void RelocEvaluation::setLowerLimits(vector<int> &lowerLimits, int &main_widx){
    //lower truncation limits are not adjusted (for now).
    
    cout << "Lower truncation cap. limits" << endl;
    for (int pidx=0; pidx<nWards; pidx++){
        cout << lowerLimits[pidx] << " " << flush;
    }
    cout << endl;
    
}

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

void RelocEvaluation::initializeStateDistribution(HeuristicQueue &hqueue){
    //initial state distribution
    
    //random number generation
    mt19937 rgen(seed);
    uniform_real_distribution<> dis(1,hqueue.Ns);
    
    pi.resize(hqueue.Ns,0);
    
    double sm=0;
    for (int i=0; i<hqueue.Ns; i++){
        pi[i] = dis(rgen);
        sm += pi[i]; 
    }
    for (int i=0; i<hqueue.Ns; i++){
        pi[i] /= sm;
    }
    
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

void RelocEvaluation::calculateWardStateSpaceSize(int ward, int numberOfWards){
    
    return((wards_pointer + ward)->calculateWardStateSpace(numberOfWards));
}

int RelocEvaluation::getWardStateSpaceSize(int ward){
    
    return((wards_pointer + ward)->wardStateSpaceSize);
}

int RelocEvaluation::getWardCapacity(int ward){
    
    return((wards_pointer + ward)->capacity);
}