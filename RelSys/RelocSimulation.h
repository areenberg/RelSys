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
#include "QueuePerformance.h"

#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;


class RelocSimulation {
public:
  
    //SIMULATION METHODS AND VARIABLES
    void setSeed(int seed);
    void simulate(double bIn, double minTime, 
        vector<int> maxWardSamples=vector<int>(1,-1),int minSamples=50); 
    void selectLogNormalServiceTime(double mult=1.0); //selects the log-normal distribution for service times 
    void disableTimeSampling(); //disables sampling of time-window sampling
    void wilsonScoreInterval(double &wilsonSpan, int &j, int &n); //calculates Wilson-score confidence interval span
    bool wilcoxonRankSum(vector<double> x, vector<double> y); //conducts the Wilcoxon rank-sum test
    void enableTimeDependency(QueuePerformance * qP); //account for time-dependent arrival rates
            
    void setAccuracy(double a); //set the tolerance of the auto simulation time 
    void setAccPref(bool ap); //set the sampling method in auto simulation time (true: all preferred arrivals, false: all arrivals to each ward)
    
    vector<vector<double>> arrivalRateMatrix; //arrival rates of each ward-patient combination
    vector<vector<double>> openTimes; //sampled open times for each ward
    vector<vector<double>> blockedTimes; //sampled blocking times for each ward
    
    //ward-patient occupancy
    vector<vector<int>> wardFreqDist; //marginal frequency distributions (all arrivals)
    vector<vector<int>> wardFreqDistPref; //marginal frequency distributions for preferred arrivals
    vector<vector<double>> wardDenDist; //marginal density distributions (all arrivals)
    vector<vector<double>> wardDenDistPref; //marginal density distributions for preferred arrivals
    vector<double> blockingProbability; //blocking probability (all arrivals)
    vector<double> blockingProbabilityPref; //blocking probability for preferred arrivals
    vector<double> expectedOccupancy; //expected server occupancy
    vector<double> expOccFraction; //expected fraction of servers occupied
    vector<int> nWardFreq; //number of samples from the marginal distributions (all arrivals)
    vector<int> nWardFreqPref; //number of samples from the marginal distribution for preferred arrivals 
    vector<vector<vector<int>>> freqDist; //frequency distribution
    vector<double> wardLoadUpperBounds;
    vector<double> wardLoadLowerBounds;
    
    //WARD INFORMATION METHODS AND VARIABLES
    int nWards; //number of wards in the system
    
    long int runtime;
    bool timeDepEnabled;
    
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
    
    int timeIndex(double cl);
    int wardSamplesToGo();
    int minTimeSamples();
    void nextServiceIdx();
    void attemptAdmission(bool &succeeded);
    bool attemptDischarge();
    double genServiceTime(int idx); //generate a random service time for the patient
    
    void openTimeTracking(int &targetWard);
    void blockedTimeTracking(int &targetWard);
    void evaluateBurnIn();
    void evalWardLoads(); //estimates load bounds for each ward
    void checkTerminate();
    bool skipWardAccuracy(int &widx);
    
    void subsetTimeSamples(int &minSamples); //randomly limits time samples to a sub-set of size minSamples 
    
    void occupancyDistTracking(int &targetWard, int &patientType);
    
    void initFreqDenDist(); //initialize the frequency/density ward-patient distributions
    void performanceMeasures(); //derives a number of performance measures    
    
    double accuracy(); //calculates the accuracy of the density distributions using wilson score intervals
    
    void printTimeSamples();
    
    long interElapsed; //intermediate elapsed time
    bool timeSamplingEnabled, checkAccuracy, checkBurnIn, burnInInit;
    vector<int> skipAccuracy; //specifies the wards to ignore when checking accuracy
    vector<vector<double>> burnInSamples;
    vector<int> bInOrder;
    vector<int> nOpenTimeSamples;
    vector<int> nBlockedTimeSamples;
    vector<vector<double>> nextArrivalTime;
    vector<double> wardStateClocks;
    vector<int> maxWrdSam;
    vector<vector<int>> wardOccupancy;
    vector<int> capUse;
    vector<double> wilsonIntervals;
    bool timeOut; //indicates if the simulation has timed out
    bool serTimeExponential; //if true, then service times are exponential; other log-normal
    double stdMult; //modifies the standard deviation in the log-normal random generator
    double accTol; //density distribution accuracy (default: 5e-3)
    bool accPref; //indicates if the preferred distribution should be used to validate accuracy (default: false) 
    double clockDis; //clock at most recent discharge;
    double wilsonSpan, wilsonMax;
    int bInSize, disIdx, cycleLen;
    
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
    //vector<Customer> nextArrival;
    Customer * nextArrival;
    QueuePerformance * qPer;
    
    //WARD INFORMATION METHODS AND VARIABLES
    QueueData * wards_pointer;
    
    //methods
    int getWardID(int ward);
    double getWardArrivalRate(int ward);
    double getWardServiceRate(int ward);
    int getWardCapacity(int ward);
    vector<double> getWardRelocationProbabilities(int ward);
    int getWardStateSpaceSize(int ward);
    double getWardTimeDep(int ward, int timeIndex);
    void calculateWardStateSpaceSize(int ward, int numberOfWards);
    void updateArrivalTime_Stationary(int widx, int pidx);
    void updateArrivalTime_TimeDep(int widx, int pidx);
    
    //LOG-NORMAL SIMULATION METHODS
    double logNormalCdfInv(double cdf,double mu,double sigma);
    double normalCdfInv(double cdf, double mu, double sigma);
    double normal01CdfInv(double p);
    double r8polyValueHorner(int m,double c[],double x);
    
};

#endif /* RELOCSIMULATION_H */