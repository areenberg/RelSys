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
#include <random>

using namespace std;


class RelocSimulation {
public:
  
    //SIMULATION METHODS AND VARIABLES
    void setSeed(int seed);
    void simulate(double bIn, double minTime, int minSamples=50); 
    
    double randomLogNormal(double mean, double std);
    
    vector<vector<double>> arrivalRateMatrix; //arrival rates of each ward-patient combination
    vector<vector<double>> openTimes; //sampled open times for each ward
    vector<vector<double>> blockedTimes; //sampled blocking times for each ward
    
    //ward-patient occupancy
    vector<vector<int>> wardFreqDist; //marginal frequency distribution
    vector<vector<vector<int>>> freqDist; //frequency distribution
    vector<vector<vector<double>>> denDist; //density distribution
    
    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system

    //CONSTRUCTOR
    //dummy constructor (not included in cpp-file) 
    RelocSimulation() {};
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

    double randomUniform(); //generate a random uniform double in the interval (0,1)
    double randomExponential(double rate); //generate a random exponentially distributed double
    //double randomLogNormal(double mean, double std);
    vector<int> randomIndices(int &from, int &to, int &len);
    
    int minTimeSamples();
    int nextServiceIdx(int &inService);
    void updateServiceArray(int idx, int &inService);
    bool attemptAdmission(int &arrIdx, vector<int> &capUse,
        vector<vector<int>> &wardOccupancy, int &inService);
    bool attemptDischarge(int &serIdx, int &inService, vector<int> &capUse);
    void updateOccupancy(vector<int> &capUse, vector<vector<int>> &wardOccupancy, int &inService);
    
    void openTimeTracking(vector<int> &nOpenTimeSamples,
 int &targetWard, vector<int> &capUse);
    
    void blockedTimeTracking(vector<int> &nBlockedTimeSamples,
 int &targetWard, vector<int> &capUse);
    
    void subsetTimeSamples(int &minSamples); //randomly limits time samples to a sub-set of size minSamples 
    
    void occupancyDistTracking(vector<vector<int>> &wardOccupancy,
    vector<int> &capUse, int &targetWard, int &patientType);
    
    void initFreqDenDist(); //initialize the frequency/density ward-patient distributions
    void freqToDensity(); //calculates the occupancy density distribution
    
    void printTimeSamples();
    
    vector<vector<double>> nextArrivalTime;
    vector<double> wardStateClocks;
    

    double burnIn, clock; //burn-in time and the simulation clock
    int simSeed; //seed for the simulation
    mt19937 rgen; //random generator
    uniform_real_distribution<> dis;

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

    //LOG-NORMAL SIMULATION METHODS
    double logNormalCdfInv(double cdf,double mu,double sigma);
    double normalCdfInv(double cdf, double mu, double sigma);
    double normal01CdfInv(double p);
    double r8polyValueHorner(int m,double c[],double x);
    
};

#endif /* RELOCSIMULATION_H */