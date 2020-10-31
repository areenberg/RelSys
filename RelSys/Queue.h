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
 * File:   Queue.h
 * Author: Anders Reenberg Andersen
 *
 * Created on September 14, 2020, 7:41 PM
 */

#ifndef QUEUE_H
#define QUEUE_H

#include "HyperQueue.h"
#include "combinatorial.h"

#include <vector>

using namespace std;

class Queue {
public:
    
    //variables
    double arrivalRate, serviceRate;
    int Ns, cap, Nh, maxNz;
    
    vector<int> state; //the current state
    int sidx; //index of the current state
    vector<int> jumpToIdx; //indices of all subsequent states (outgoing transition)
    vector<int> jumpFromIdx; //indices of all previous states (ingoing transition)
    vector<double> jumpToRate, jumpFromRate;
    int toIdxSize, fromIdxSize; //number of values allocated to jumpToIdx and jumpFromIdx, respectively

    vector<vector<int>> qColumnIndices;
    vector<vector<double>> qValues;

    //methods
    void buildChain();
    void buildTransposedChain();
    void printChain();
    
    void initializeState(); //initialize the current state
    void nextCurrentState(); //move to the next state
    
    void allOutgoing(); //generate all jumpToIdx relative to the current state      
    void allIngoing(); //generate all jumpFromIdx relative to the current state
    
    vector<double> marginalDist(vector<double> &pi); //returns the marginal distribution
    double expectedOccupancy(vector<double> &pi); //returns the expected occupancy of capacity
    double expectedOccupancyFraction(vector<double> &pi); //returns the expected fraction of capacity utilized
    double rejectionProbability(vector<double> &pi); //returns the rejection probability using the marginal occupancy distribution

    //constructor and destructor
    Queue(int c, vector<int> upperLim, vector<int> lowerLim, double aRate, double sRate, int nhq, HyperQueue* hbQueues);
    Queue(const Queue& orig);
    virtual ~Queue();

private:

    //parameters and pointers
    HyperQueue * hbQueues_pointer;
    
    vector<int> upperLim; //capacity limits
    vector<int> lowerLim;
    
    combinatorial cmb; //used for calculating jumps
    
    int K_use;

    //methods
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
    
};

#endif /* QUEUE_H */

