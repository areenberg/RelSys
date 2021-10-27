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
 * File:   RelocEvaluation.h
 * Author: Anders Reenberg Andersen
 *
 * Created on November 24, 2020, 9:31 PM
 */

#ifndef RELOCEVALUATION_H
#define RELOCEVALUATION_H

#include "HeuristicQueue.h"
#include "HyperQueue.h"
#include "QueueData.h"
#include "RelocSimulation.h"


#include <vector>

using namespace std;

class RelocEvaluation {
public:
    
    //METHODS
    void runHeuristic(int main_widx);
    void runSimulation(int seed, int burnIn,
        int minTime, int minSamples=50);
    void simulateMarginalDist(double sampleBurnIn, int collectSamples); //replace numerical evaluation with a simulation
    void setBinMap(vector<vector<int>> bMap);
    void setBinDischargeRates(vector<double> disRates);
    void setOpenHyperStates(int n);
    void setBlockedHyperStates(int n);
    
    //VARIABLES
    vector<double> pi; //state distribution
    vector<double> marginalDist; //marginal distribution of the main queue
    double blockingProbability, expectedOccupancy, expOccFraction;
    double vmemory; //virtual memory consumption during runtime in kilobytes
    
    RelocEvaluation() {}; //dummy constructor 
    RelocEvaluation(int nW, QueueData * wards);
    RelocEvaluation(const RelocEvaluation& orig);
    virtual ~RelocEvaluation();
private:

    //METHODS
    void initializeSystem();
    void setDefaultBinMap();
    void setDefaultHyperQueueStates();
    void initializeStateDistribution(HeuristicQueue &hqueue, bool erlangInit=false);
    void setUpperLimits(vector<int> &upperLimits, int &main_widx);
    void setLowerLimits(vector<int> &lowerLimits, int &main_widx);
    
    double chebyshevBound(double &x, vector<int> &freqDist);
    double Gfunction(double q, double x);
    double sampleMean(vector<int> &freqDist);
    double sampleSD(vector<int> &freqDist);
    double erlangLoss(int &k, double &lambda, double &mu, int &servers);
    long double factorial(int &x);
    
    //VARIABLES
    vector<vector<int>> binMap;
    vector<double> newDisRates;
    vector<int> hyperOpenStates; //number of open states for each hyper queue
    vector<int> hyperBlockedStates; //number of blocked states for each hyper queue
    bool simReady,simMargDist, changeDisRates;
    int nWards, seed, collectSamples;
    double sampleBurnIn;
    RelocSimulation * sim_pointer;
    
    mt19937 rgen; //random generator
    uniform_real_distribution<> dis;
    
    
    //WARD INFORMATION METHODS AND VARIABLES
    QueueData * wards_pointer;
    
    //methods
    int getWardID(int ward);
    double getWardArrivalRate(int ward);
    double getWardServiceRate(int ward);
    int getWardCapacity(int ward);
    vector<double> getWardRelocationProbabilities(int ward);
    int getWardStateSpaceSize(int ward);
    void calculateWardStateSpaceSize(int ward, int numberOfWards);

    
};

#endif /* RELOCEVALUATION_H */

