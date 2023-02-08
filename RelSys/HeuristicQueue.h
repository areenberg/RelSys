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
 * File:   HeuristicQueue.h
 * Author: Anders Reenberg Andersen
 *
 * Created on September 14, 2020, 7:41 PM
 */

#ifndef HEURISTICQUEUE_H
#define HEURISTICQUEUE_H

#include "HyperQueue.h"
#include "Combinatorial.h"

#include <vector>

using namespace std;

class HeuristicQueue {
public:
    
    //variables
    double arrivalRate,serviceRate;
    int Ns,cap,Nh,maxNz,nBins;
    
    vector<int> state; //the current state
    int sidx; //index of the current state
    vector<int> jumpToIdx; //indices of all subsequent states (outgoing transition)
    vector<int> jumpFromIdx; //indices of all previous states (ingoing transition)
    vector<double> jumpToRate, jumpFromRate;
    int toIdxSize, fromIdxSize; //number of values allocated to jumpToIdx and jumpFromIdx, respectively

    vector<vector<int>> qColumnIndices;
    vector<vector<double>> qValues;
    vector<double> margDist; //marginal probability distribution

    //methods
    void buildChain();
    void buildTransposedChain();
    void printStateSpace();
    void newbinDischargeRates(vector<double> newBinDisRates);
    
    void initializeState(); //initialize the current state
    void nextCurrentState(); //move to the next state
    void previousCurrentState(); //move to the previous state
    
    void allOutgoing(); //generate all jumpToIdx relative to the current state      
//    void allIngoing(); //generate all jumpFromIdx relative to the current state (currently buggy)
    
    void marginalDist(vector<double> &pi); //returns the marginal distribution
    double expectedOccupancy(); //returns the expected occupancy of capacity
    double rejectionProbability(); //returns the rejection probability using the marginal occupancy distribution

    //constructor and destructor
    HeuristicQueue() {}; //dummy constructor (not included in cpp-file)
    HeuristicQueue(int main_widx, vector<vector<int>> binMap, int c, vector<int> upperLim, vector<int> lowerLim,
    double aRate, double sRate, int nhq, HyperQueue* hbQueues, QueueData * wards);
    HeuristicQueue(const HeuristicQueue& orig);
    virtual ~HeuristicQueue();

private:

    //parameters and pointers
    HyperQueue * hbQueues_pointer;
    QueueData * wards_pointer;
    
    vector<vector<int>> binMap; //patient types x bin map
    vector<int> upperLim; //capacity limits
    vector<int> lowerLim;
    vector<double> binDischargeRates;
    vector<int> newHypIdx;
    vector<int> newWardIdx;
    int main_widx; //index of the main ward
    int maximumSize;
    vector<int> hyperWidx_vector;
    vector<int> fwvec;
    vector<bool> hblocked;
    
    Combinatorial cmb; //used for calculating jumps
    
    int K_use;

    //methods
    double getAdmissionRate(int &bidx, vector<bool> &hblocked, bool &hasMain);
    void calculateDischargeRates();
    vector<int> redundantWards(); //find the wards that are redundant to the main ward
    void optimizeRelocNetwork();
    void adjustLimits();
    void calculateSize(); //calculate and store size of the state space
    void initializeJumbVectors(); //allocate memory for jumpToIdx and jumpFromIdx
    void maximumNonZero(); //derives an upper bound for the number of non-zero elements in the transition rate matrix;
    double calculateDiagonal(); //calculates the diagonal related to the current state
    void checkInput();
    int forwardOne(int Ku, vector<int> &j, int targetval, int targetidx);
    int backwardOne(int Ku, vector<int> &j, int targetval, int targetidx);


    int hyperOpenStates(int hq); //get number of open states in hyperqueue hq
    int hyperBlockedStates(int hq); //get number of blocked states in hyperqueue hq
    double getHyperOpenRate(int hq, int idx);
    double getHyperOpenDist(int hq, int idx);
    double getHyperBlockedRate(int hq, int idx);
    double getHyperBlockedDist(int hq, int idx);
    int getHyperSize(int hq);
    double getHyperArrivalRate(int hq);
    double getHyperServiceRate(int hq);
    double getWardArrivalRate(int widx);
    double getWardServiceRate(int widx);
    vector<double> getWardRelocationProbabilities(int widx);
    
};

#endif /* HEURISTICQUEUE_H */

