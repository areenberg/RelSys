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
 * File:   RelocSimulation.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on November 13, 2020, 2:02 PM
 */

#include "RelocSimulation.h"
#include "Patient.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <limits>
#include <random>

using namespace std;

RelocSimulation::RelocSimulation(int nW, WardData * wards):
wards_pointer(wards),
nWards(nW)        
{
    initializeSystem();
}

RelocSimulation::RelocSimulation(const RelocSimulation& orig) {
}

RelocSimulation::~RelocSimulation() {
}

void RelocSimulation::initializeSystem(){
    calculateArrivalRates();
    
    //setup for the random number generator
    random_device rd;
    simSeed = rd();
    mt19937 rgen(simSeed);
    uniform_real_distribution<> dis(0,1);
}

void RelocSimulation::calculateArrivalRates(){
    //calculates the arrival rate into each bin of the system
    
    arrivalRateMatrix.resize(nWards);
    for (int widx=0; widx<nWards; widx++){
        arrivalRateMatrix[widx].resize(nWards,0);
        for (int pidx=0; pidx<nWards; pidx++){
            if (widx==pidx){
                arrivalRateMatrix[widx][pidx] = getWardArrivalRate(widx);
            }else{
                arrivalRateMatrix[widx][pidx] = getWardArrivalRate(pidx)*
                        getWardRelocationProbabilities(pidx)[widx];
            }
        }
    }
    
}

void RelocSimulation::setSeed(int seed){
    //srand(seed);
    simSeed = seed;
    mt19937 rgen(simSeed);
}

void RelocSimulation::initializeArrivalTimes(double currentClock){
    //allocate memory and fill in the entire nextArrivalTime matrix
    
    nextArrivalTime.clear();
    nextArrivalTime.resize(nWards);
    for (int i=0; i<nWards; i++){
        nextArrivalTime[i].resize(nWards,currentClock);
    }
    
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=0; pidx<nWards; pidx++){
            nextArrivalTime[widx][pidx] += randomExponential(arrivalRateMatrix[widx][pidx]);
        }
    }
    
}

void RelocSimulation::updateArrivalTime(int widx, int pidx){
    //derives the next arrival time for a single ward-patient pair
    
    nextArrivalTime[widx][pidx] += randomExponential(arrivalRateMatrix[widx][pidx]);
    
}

void RelocSimulation::generateArrivalList(int length, double currentClock){
    //generate a long list of future arrivals with arrival times offset
    //by currentClock
    
    initializeArrivalTimes(currentClock);
    
    //append patients to the list
    double arrClock, serTime, mn;
    int min_widx, min_pidx;
    //vector<int> genWidx(nWards,0);
    //vector<int> genPidx(nWards,0);
    
    for (int i=0; i<length; i++){
        
        mn = numeric_limits<double>::max();
        
        //find minimum arrival time        
        for (int widx=0; widx<nWards; widx++){
            for (int pidx=0; pidx<nWards; pidx++){
                if (nextArrivalTime[widx][pidx]<mn){
                    mn = nextArrivalTime[widx][pidx];
                    min_widx = widx; min_pidx = pidx;
                }
            }
        }        
        
        
        //append patient
        arrClock = nextArrivalTime[min_widx][min_pidx];
        serTime = randomExponential(getWardServiceRate(min_pidx));
//        serTime = 0.0;
        arrival_array[i] = Patient(arrClock,serTime,min_widx,min_pidx);
        
        //update ward-patient
        updateArrivalTime(min_widx,min_pidx);
    }
    
    
}

void RelocSimulation::simulate(double bIn, double minTime, int minSamples){
    //burnIn is the burn-in time.
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
    //variables
    double nextArrClock, nextSerClock;
    int patientArraySize, maxOcc, inService, arrIdx, 
            serIdx, targetWard;
    bool succeeded, note = false;
    
    //vector<vector<int>> wardOccupancy;
    vector<int> nOpenTimeSamples(nWards,0);
    vector<int> nBlockedTimeSamples(nWards,0);
    wardStateClocks.resize(nWards,0);
    
    //--------------------
    //Initialize
    //--------------------
    
    //burn-in time
    burnIn = bIn;
    
    //simulation clock
    clock = 0;
    
    //patient arrivals
    patientArraySize = 1e5;
    arrival_array = new Patient[patientArraySize];
    generateArrivalList(patientArraySize,clock);
    arrIdx = 0;
    
    //patients in service
    maxOcc = 0; inService = 0;
    for (int widx=0; widx<nWards; widx++){
        maxOcc += getWardCapacity(widx);
    }
    service_array = new Patient[maxOcc];
    
    //occupancy of the system
    vector<int> capUse(nWards,0);
    
    //tracking
    openTimes.resize(nWards);
    blockedTimes.resize(nWards);
    
    //-------------------
    //Simulate
    //-------------------
    cout << "Running simulation..." << flush;
    while (arrIdx<patientArraySize && (clock<minTime || 
            minTimeSamples()<minSamples)){
        
        //clock at next arrival
        nextArrClock = arrival_array[arrIdx].arrivalClock;
        
//        cout << "sim. clock = " << clock << endl;
//        for (int i=0; i<nWards; i++){
//            cout << wardStateClocks[i] << ", ";
//        }
//        cout << endl;
//        cout << "admission widx=" << arrival_array[arrIdx].wardTarget << ", pidx=" << arrival_array[arrIdx].patientType << ", adm.clock=" << nextArrClock << endl; 
//        cout << "cap use=" << capUse[arrival_array[arrIdx].wardTarget] << ", cap lim=" << getWardCapacity(arrival_array[arrIdx].wardTarget) << endl;
//    
//        cout << "old array:" << endl;
//        for (int i=0; i<inService; i++){
//            cout << service_array[i].serviceClock << "," << service_array[i].wardTarget << "," << service_array[i].patientType << endl;    
//        }
    
        
        if (inService>0){
            //clock at next discharge
            serIdx = nextServiceIdx(inService);
            nextSerClock = service_array[serIdx].serviceClock;
//            cout << "inService=" << inService << ", serIdx=" << serIdx << ", nextSerClock=" << nextSerClock << endl;
            
            if (nextSerClock<nextArrClock){
                clock = nextSerClock;
                targetWard = service_array[serIdx].wardTarget;
                succeeded = attemptDischarge(serIdx,inService,capUse);
                if (succeeded){
                    blockedTimeTracking(nBlockedTimeSamples,targetWard,capUse);
                }
            }else{
                clock = nextArrClock;
                succeeded = attemptAdmission(arrIdx,capUse,inService);
                if (succeeded){
                    openTimeTracking(nOpenTimeSamples,
                        service_array[inService-1].wardTarget,capUse);
                }
            }
        }else{
            clock = nextArrClock; 
            succeeded = attemptAdmission(arrIdx,capUse,inService);
            if (succeeded){
                openTimeTracking(nOpenTimeSamples,
                    service_array[inService-1].wardTarget,capUse);
            }
        }
        
        updateOccupancy(capUse,inService);
    
        
//        cout << "new array:" << endl;
//        for (int i=0; i<inService; i++){
//            cout << service_array[i].serviceClock << "," << service_array[i].wardTarget << "," << service_array[i].patientType << endl; 
//        }
//        cout << "----------" << endl;
    
        
        //generate a new batch of patient arrivals
        if (arrIdx==(patientArraySize-1)){
            generateArrivalList(patientArraySize,clock);
            arrIdx = 0;
        }else if(note==false && clock>minTime && minTimeSamples()<minSamples){
            cout << "still collecting samples..." << flush;
            note = true;
        }else if(note==false && clock<minTime && minTimeSamples()>=minSamples){
            cout << "state-time samples collected..." << flush;
            note = true;
        }
        
    }
    cout << "done." << endl;
    
    printTimeSamples();
    
}

void RelocSimulation::printTimeSamples(){
    
    cout << "\n----------------------------" << endl;
    cout << "*Number of open-time samples*" << endl;
    for (int i=0; i<openTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << openTimes[i].size();
        if (i<(openTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << "\n" << endl;
    cout << "*Number of blocked-time samples*" << endl;
    for (int i=0; i<blockedTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << blockedTimes[i].size();
        if (i<(blockedTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << endl;
    cout << "----------------------------\n" << endl;
    
}

void RelocSimulation::updateOccupancy(vector<int> &capUse, int &inService){
    
    capUse.clear();
    capUse.resize(nWards,0);
    
    for (int i=0; i<inService; i++){
        capUse[service_array[i].wardTarget]++;
    }
    
}

void RelocSimulation::openTimeTracking(
        vector<int> &nOpenTimeSamples, int &targetWard, vector<int> &capUse){
    
    if ((capUse[targetWard]+1)==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
            openTimes[targetWard].push_back(sample);
            nOpenTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}
    
void RelocSimulation::blockedTimeTracking(
        vector<int> &nBlockedTimeSamples, int &targetWard, vector<int> &capUse){
    
    if (capUse[targetWard]==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
//            if (sample==0){
//                cout << "SAMPLE 0" << endl;
//            }
            blockedTimes[targetWard].push_back(sample);
            nBlockedTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}


bool RelocSimulation::attemptDischarge(int &serIdx, int &inService, vector<int> &capUse){
    
    //subtract patient from service
    inService--;
    
    updateServiceArray(serIdx,inService);
    
    //cout << "attempted discharge success" << endl;
    
    return(true);
}

bool RelocSimulation::attemptAdmission(int &arrIdx, vector<int> &capUse, int &inService){
    
    bool Ok = false;
    if (arrival_array[arrIdx].wardTarget==arrival_array[arrIdx].patientType &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }else if(capUse[arrival_array[arrIdx].patientType]==getWardCapacity(arrival_array[arrIdx].patientType) &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }
    
    if (Ok){
        
        //insert into service array
        service_array[inService] = Patient(arrival_array[arrIdx].arrivalClock,arrival_array[arrIdx].serviceTime,
                arrival_array[arrIdx].wardTarget,arrival_array[arrIdx].patientType);
        //calculate and insert clock at discharge
        service_array[inService].serviceClock = clock + service_array[inService].serviceTime;
        
        //adjust number of patients currently in service
        inService++;
    }
    
    //move to next arrival
    arrIdx++;
    
    return(Ok);
}
    
void RelocSimulation::updateServiceArray(int idx, int &inService){
    //removes a patient in service and adjusts the entire list
    //assumes the patient <<has been subtracted>> from inService
    
    if (idx<inService){
        double arrClock, serTime, serClock;
        int widx, pidx;
        
        for (int i=idx; i<inService; i++){
            arrClock = service_array[i+1].arrivalClock;
            serTime = service_array[i+1].serviceTime;
            serClock = service_array[i+1].serviceClock;
            widx = service_array[i+1].wardTarget;
            pidx = service_array[i+1].patientType;
            
            service_array[i] = Patient(arrClock,serTime,widx,pidx);
            service_array[i].serviceClock = serClock;
        }
    }
    
}

int RelocSimulation::nextServiceIdx(int &inService){
    
    if (inService>1){
        int idx;
        double mn = numeric_limits<double>::max();
        for (int i=0; i<inService; i++){
            if (service_array[i].serviceClock<mn){
               mn = service_array[i].serviceClock;
               idx = i;
            }
        }
//        cout << "serIdx=" << idx << ", clock=" << mn << endl;
        return(idx);
    }else{
        return(0);
    }
    
}

int RelocSimulation::minTimeSamples(){
    
    int mn = numeric_limits<int>::max();
    
    for (int i=0; i<openTimes.size(); i++){
        if (mn>openTimes[i].size()){
            mn = openTimes[i].size();
        }
    }
    for (int i=0; i<blockedTimes.size(); i++){
        if (mn>blockedTimes[i].size()){
            mn = blockedTimes[i].size();
        }
    }
    
    return(mn);
}


double RelocSimulation::randomExponential(double rate){
    //generate a random uniform double in the interval (0,1)
    return(log(1-randomUniform())/(-rate));
}

double RelocSimulation::randomUniform(){
    //double r = (double)rand() / RAND_MAX;
    double r = dis(rgen);
    return(r);
}

int RelocSimulation::getWardID(int ward){
    
    return((wards_pointer + ward)->wardnr);
}

double RelocSimulation::getWardArrivalRate(int ward){
    
    return((wards_pointer + ward)->arrivalRate);
}

double RelocSimulation::getWardServiceRate(int ward){
    
    return((wards_pointer + ward)->serviceRate);
}

vector<double> RelocSimulation::getWardRelocationProbabilities(int ward){
    
    return((wards_pointer + ward)->relocationProbabilities);
}

void RelocSimulation::calculateWardStateSpaceSize(int ward, int numberOfWards){
    
    return((wards_pointer + ward)->calculateWardStateSpace(numberOfWards));
}

int RelocSimulation::getWardStateSpaceSize(int ward){
    
    return((wards_pointer + ward)->wardStateSpaceSize);
}

int RelocSimulation::getWardCapacity(int ward){
    
    return((wards_pointer + ward)->capacity);
}