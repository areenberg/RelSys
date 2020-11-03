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

#include "WardData.h"

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
    void calculateStateSpaceSize();
    void buildTransposedChain();
    void initializeState(); //initialize the current state
    void nextCurrentState(); //move to the next state
    void printState();
    
    void allIngoing(); //generate all jump indices and rates relative to the current state 

    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system

    //CONSTRUCTOR
    EntireSystem(int nW, WardData * wards);
    EntireSystem(const EntireSystem& orig);
    virtual ~EntireSystem();

private:

    //MARKOV CHAIN METHODS
    void initializeSystem();
    void initializeJumbVectors();

    //WARD INFORMATION METHODS AND VARIABLES
    WardData * wards_pointer;
    
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