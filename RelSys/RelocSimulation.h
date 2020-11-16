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
 * File:   RelocSimulation.h
 * Author: Anders Reenberg Andersen
 * 
 */

#ifndef RELOCSIMULATION_H
#define RELOCSIMULATION_H

#include "WardData.h"
#include "Patient.h"

#include <iostream>
#include <vector>

using namespace std;


class RelocSimulation {
public:
  
    //SIMULATION METHODS AND VARIABLES
    void setSeed(int seed);
    void simulate(double burnIn, double minTime, int minSamples=50); 

    vector<vector<double>> arrivalRateMatrix; //arrival rates of each ward-patient combinatiom
    
    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system

    //CONSTRUCTOR
    RelocSimulation(int nW, WardData * wards);
    RelocSimulation(const RelocSimulation& orig);
    virtual ~RelocSimulation();

private:

    //SIMULATION METHODS AND VARIABLES
    void initializeSystem();
    void calculateArrivalRates();
    void generateArrivalList(int length, double currentClock);
    void updateArrivalTime(int widx, int pidx); //update the nextArrivalTime matrix
    void initializeArrivalTimes(double currentClock); //initialize the entire nextArrivalTime matrix

    double randomUniform(double from, double to); //generate a random double between from and to
    double randomExponential(double rate); //generate a random exponentially distributed double

    int minTimeSamples(vector<int> &openTimes, vector<int> &blockedTimes);
    int nextServiceIdx(int &inService);
    void updateServiceArray(int idx, int &inService);
    bool attemptAdmission(double &currentClock, int &arrIdx, vector<vector<int>> &wardOccupancy, vector<int> &capUse, int &inService);
    bool attemptDischarge(double & currentClock, int &serIdx, int &inService, vector<vector<int>> &wardOccupancy, vector<int> &capUse);
    
    void openTimeTracking(bool &success, double &currentClock, vector<vector<double>> &openTimes,
    vector<double> &wardStateTimes, int &targetWard, vector<int> &capUse);
    
    void blockedTimeTracking(bool &success, double &currentClock, vector<vector<double>> &blockedTimes,
    vector<double> &wardStateTimes, int &targetWard, vector<int> &capUse);
    
    vector<vector<double>> openTimes; //sampled open times for each ward
    vector<vector<double>> blockedTimes; //sampled blocking times for each ward
    
    vector<vector<double>> nextArrivalTime;
    
    //PATIENT METHODS AND VARIABLES
    Patient * arrival_array;
    Patient * service_array;
    
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

#endif /* RELOCSIMULATION_H */