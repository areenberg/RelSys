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

#include "QueueData.h"
#include "Customer.h"

#include <iostream>
#include <vector>
#include <random>

using namespace std;


class RelocSimulation {
public:
  
    //SIMULATION METHODS AND VARIABLES
    void setSeed(int seed);
    void simulate(double bIn, double minTime, 
        vector<int> maxWardSamples=vector<int>(1,-1),int minSamples=50); 
    void selectLogNormalServiceTime(double mult=1.0); //selects the log-normal distribution for service times 
    void disableTimeSampling(); //disables sampling of time-window sampling
    vector<double> wilsonScoreInterval(double p, int n); //calculates Wilson-score confidence intervals
    bool wilcoxonRankSum(vector<double> x, vector<double> y); //conducts the Wilcoxon rank-sum test
    
    void setAccuracy(double a);
    
    
    vector<vector<double>> arrivalRateMatrix; //arrival rates of each ward-patient combination
    vector<vector<double>> openTimes; //sampled open times for each ward
    vector<vector<double>> blockedTimes; //sampled blocking times for each ward
    
    //ward-patient occupancy
    vector<vector<int>> wardFreqDist; //marginal frequency distributions
    vector<vector<double>> wardDenDist; //marginal density distributions
    vector<double> blockingProbability; //blocking probability
    vector<double> expectedOccupancy; //expected server occupancy
    vector<double> expOccFraction; //expected fraction of servers occupied
    vector<int> nWardFreq; //number of samples in the marginal distributions
    vector<vector<vector<int>>> freqDist; //frequency distribution
    vector<vector<vector<double>>> denDist; //density distribution
    
    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system
    
    long int runtime;
    
    //CONSTRUCTOR
    //dummy constructor (not included in cpp-file) 
    RelocSimulation() {};
    RelocSimulation(int nW, QueueData * wards);
    RelocSimulation(const RelocSimulation& orig);
    virtual ~RelocSimulation();

private:

    //SIMULATION METHODS AND VARIABLES
    void initializeSystem();
    void calculateArrivalRates();
    void generateArrivalList(double currentClock);
    void generateArrival();
    void updateArrivalTime(int widx, int pidx); //update the nextArrivalTime matrix
    void initializeArrivalTimes(double currentClock); //initialize the entire nextArrivalTime matrix

    double randomUniform(); //generate a random uniform double in the interval (0,1)
    double randomExponential(double rate); //generate a random exponentially distributed double
    double randomLogNormal(double mean, double std);
    vector<int> randomIndices(int &from, int &to, int &len);
    
    int wardSamplesToGo();
    int minTimeSamples();
    void nextServiceIdx();
    bool attemptAdmission(int &arrIdx);
    bool attemptDischarge();
    void updateOccupancy();
    double genServiceTime(int idx); //generate a random service time for the patient
    
    void openTimeTracking(int &targetWard);
    void blockedTimeTracking(int &targetWard);
    void evaluateBurnIn();
    
    void subsetTimeSamples(int &minSamples); //randomly limits time samples to a sub-set of size minSamples 
    
    void occupancyDistTracking(int &targetWard, int &patientType);
    
    void initFreqDenDist(); //initialize the frequency/density ward-patient distributions
    void freqToDensity(); //calculates the occupancy density distribution
    void performanceMeasures(); //derives a number of performance measures    
    
    double accuracy(); //calculates the accuracy of the density distributions using wilson score intervals
    
    void printTimeSamples();
    
    bool timeSamplingEnabled, checkAccuracy, checkBurnIn, burnInInit;
    vector<vector<double>> burnInSamples;
    vector<int> bInOrder;
    vector<int> nOpenTimeSamples;
    vector<int> nBlockedTimeSamples;
    vector<vector<double>> nextArrivalTime;
    vector<double> wardStateClocks;
    vector<int> maxWrdSam;
    vector<vector<int>> wardOccupancy;
    vector<int> capUse;
    bool serTimeExponential; //if true, then service times are exponential; other log-normal
    double stdMult; //modifies the standard deviation in the log-normal random generator
    double accTol; //density distribution accuracy
    double clockDis; //clock at most recent discharge;
    int bInSize, disIdx;
    
    int min_widx, min_pidx;
    double simTime, burnIn, clock; //sim. time, burn-in time and the simulation clock
    int simSeed; //seed for the simulation
    int inService; //number of customers in service
    int serIdx,insIdx; //indices for next service and recently inserted
    int patientArraySize, maxOcc; //size of patient arrival and service arrays
    //mt19937 rgen; //random generator
    //default_random_engine rgen;
    default_random_engine rgen;
    uniform_real_distribution<> dis;
    int mxRnd;
    
    //PATIENT METHODS AND VARIABLES
    vector<Customer> service_array;
    vector<Customer> nextArrival;
    
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

    //LOG-NORMAL SIMULATION METHODS
    double logNormalCdfInv(double cdf,double mu,double sigma);
    double normalCdfInv(double cdf, double mu, double sigma);
    double normal01CdfInv(double p);
    double r8polyValueHorner(int m,double c[],double x);
    
};

#endif /* RELOCSIMULATION_H */