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
#include "WardData.h"
#include "RelocSimulation.h"
#include "PhaseFitter.h"

#include <vector>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

RelocEvaluation::RelocEvaluation(int nW, WardData * wards):
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
    
    sim_pointer[0].simulate(burnIn,minTime,minSamples);
    
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
        HyperQueue * hq_array = new HyperQueue[nhq];
        for (int i=0; i<nhq; i++){
            arr =  getWardArrivalRate(hyperWidx_vector[i])*getWardRelocationProbabilities(hyperWidx_vector[i])[main_widx];
            
            hq_array[i] = HyperQueue(hyperWidx_vector[i],statesBlocked,statesOpen,
                arr,getWardServiceRate(hyperWidx_vector[i]),sim_pointer);
            
            hq_array[i].fitAll(seed);
        }
        
        
        //specifications of main queue
        double arrivalRate = getWardArrivalRate(main_widx);
        double serviceRate = getWardServiceRate(main_widx);
        int capacity = getWardCapacity(main_widx);
    
        //upper and lower limits in each queue. must include all queues (main and hyper queues)
        vector<int> upperLimits(nWards,0);
        vector<int> lowerLimits(nWards,0);
        //automatically adjust truncation
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
        
        //store the marginal distribution and some other metrics
        marginalDist = hqueue.marginalDist(pi);
        expectedOccupancy = hqueue.expectedOccupancy(pi);
        expOccFraction = hqueue.expectedOccupancyFraction(pi);
        blockingProbability = marginalDist[marginalDist.size()-1];
        
        
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        cout << "Runtime of CTMC (excl. simulation): " << duration.count() << " milliseconds\n" << endl;
        
        //print and store some general metrics
        cout << "\nWARD " << (main_widx+1) << " STATS:" << endl;
        
        cout << "Expected load = " << expectedOccupancy << endl;
        cout << "Capacity utilization = " << expOccFraction*100 << "%\n" << endl;
        
        
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
    double tailden = 1e-4;
    
    vector<vector<double>> dd = sim_pointer[0].denDist[main_widx];
    double sm;
    int k;
    
    cout << "Upper truncation cap. limits" << endl;
    for (int pidx=0; pidx<nWards; pidx++){
        k = getWardCapacity(main_widx); sm = dd[pidx][k];
        while (k>0 && sm<tailden){
            k--;
            sm += dd[pidx][k];
        }
        upperLimits[pidx] = k;
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