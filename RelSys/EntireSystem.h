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
 * File:   EntireSystem.h
 * Author: Anders Reenberg Andersen
 * 
 */

#ifndef ENTIRESYSTEM_H
#define ENTIRESYSTEM_H

#include "QueueData.h"

#include <iostream>
#include <vector>

using namespace std;

class EntireSystem {
public:
    
    //MARKOV CHAIN METHODS AND VARIABLES
    vector<vector<int>> state; //the current state (wards x patients)
    int sidx; //index of the current state
    int nS; //state space size for the entire system
    vector<int> jumpFromIdx; //indices of all previous states (ingoing transition)
    vector<double> jumpFromRate;
    int fromIdxSize; //number of values allocated to jumpFromIdx

    vector<vector<int>> qColumnIndices;
    vector<vector<double>> qValues;
    
    vector<int> K_use; //current capacity use in each ward

    //methods
    void selfCheck(); //runs a self-check targeting the structure of the transition matrix
    void calculateStateSpaceSize();
    void buildTransposedChain();
    void initializeState(); //initialize the current state
    void nextCurrentState(); //move to the next state
    void printState();
    void printStateSpace();
    
    void allIngoing(); //generate all jump indices and rates relative to the current state 
    vector<double> marginalDist(int widx, vector<double> &pi);

    double expectedOccupancy(int widx, vector<double> &pi); //returns the expected occupancy of capacity
    double expectedOccupancyFraction(int widx, vector<double> &pi); //returns the expected fraction of capacity utilized
    double rejectionProbability(int widx, vector<double> &pi); //returns the rejection probability using the marginal occupancy distribution

    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system

    //CONSTRUCTOR
    EntireSystem(int nW, QueueData * wards);
    EntireSystem(const EntireSystem& orig);
    virtual ~EntireSystem();

private:

    //MARKOV CHAIN METHODS
    void initializeSystem();
    void initializeJumbVectors();
    int forwardOne(int wardCapacityUsed, vector<int> &j, int targetval, int &pidx, int &widx);
    int backwardOne(int wardCapacityUsed, vector<int> &j, int targetval, int &pidx, int &widx);
    double diagonalRate(); //value of diagonal related to the current state
    
    //self-check methods
    void nextStateLocalInput(vector<vector<int>> &s, vector<int> &capuse);
    bool checkOneStateChange(int currentStateIdx, vector<vector<int>> &currentState,vector<vector<int>> &fromState, double rate, int idx);
    void transposeChain();
    bool checkDiagonalRates();
    vector<vector<int>> qColumnIndices_t; //transposed transition matrix
    vector<vector<double>> qValues_t;
    

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

#endif /* ENTIRESYSTEM_H */