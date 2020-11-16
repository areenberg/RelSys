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
    srand(seed);
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
    arrival_array = new Patient[length];
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
        arrival_array[i] = Patient(arrClock,serTime,min_widx,min_pidx);
        
        //update ward-patient
        updateArrivalTime(min_widx,min_pidx);
    }
    
    
}

void RelocSimulation::simulate(double burnIn, double minTime, int minSamples){
    //burnIn is the burn-in time.
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
    //variables
    double clock, nextArrClock, nextSerClock;
    int patientArraySize, maxOcc, inService, arrIdx, 
            serIdx, targetWard;
    bool succeeded;
    
    vector<vector<int>> wardOccupancy;
    vector<int> nOpenTimeSamples(nWards,0);
    vector<int> nBlockedTimeSamples(nWards,0);
    vector<double> wardStateTimes(nWards,0);
    
    //--------------------
    //Initialize
    //--------------------
    
    //simulation clock
    clock = 0;
    
    //patient arrivals
    patientArraySize = 1e4;
    generateArrivalList(patientArraySize,clock);
    arrIdx = 0;
    
    //patients in service
    maxOcc = 0; inService = 0;
    for (int widx=0; widx<nWards; widx++){
        maxOcc += getWardCapacity(widx);
    }
    service_array = new Patient[maxOcc];
    
    //occupancy of the system
    wardOccupancy.resize(nWards);
    for (int widx=0; widx<nWards; widx++){
        wardOccupancy[widx].resize(nWards,0);
    }
    vector<int> capUse(nWards,0);
    
    //tracking
    openTimes.resize(nWards);
    blockedTimes.resize(nWards);
    
    //-------------------
    //Simulate
    //-------------------
    while (arrIdx<patientArraySize && (clock<minTime || 
            minTimeSamples(nOpenTimeSamples,nBlockedTimeSamples)<minSamples)){
        
        //clock at next arrival
        nextArrClock = arrival_array[arrIdx].arrivalClock;
        
        if (inService>0){
            //clock at next discharge
            serIdx = nextServiceIdx(inService);
            nextSerClock = service_array[serIdx].serviceClock;
            
            if (nextSerClock<nextArrClock){
                targetWard = service_array[serIdx].wardTarget;
                succeeded = attemptDischarge(clock,serIdx,inService,wardOccupancy,capUse);
                if (clock>burnIn){
                    blockedTimeTracking(succeeded,clock,blockedTimes,
                        wardStateTimes,targetWard,capUse);
                }
            }else{
                succeeded = attemptAdmission(clock,arrIdx,wardOccupancy,capUse,inService);
                if (clock>burnIn){
                    openTimeTracking(succeeded,clock,openTimes,wardStateTimes,
                        service_array[inService-1].wardTarget,capUse);
                }
            }
        }else{
            succeeded = attemptAdmission(clock,arrIdx,wardOccupancy,capUse,inService);
            if (clock>burnIn){
                openTimeTracking(succeeded,clock,openTimes,wardStateTimes,
                    service_array[inService-1].wardTarget,capUse);
            }
        }
        
        cout << "Sim. clock = " << clock << endl;
    }
    
}

void RelocSimulation::openTimeTracking(bool &success, double &currentClock, vector<vector<double>> &openTimes,
        vector<double> &wardStateTimes, int &targetWard, vector<int> &capUse){
    
    if (success && capUse[targetWard]==getWardCapacity(targetWard)){
        double sample = currentClock - wardStateTimes[targetWard];
        openTimes[targetWard].push_back(sample);
        wardStateTimes[targetWard] = currentClock;
        cout << "added open sample = " << sample << " to ward " << targetWard << endl;
    }
    
}
    
void RelocSimulation::blockedTimeTracking(bool &success, double &currentClock, vector<vector<double>> &blockedTimes,
        vector<double> &wardStateTimes, int &targetWard, vector<int> &capUse){
    
    if (success && (capUse[targetWard]+1)==getWardCapacity(targetWard)){
        double sample = currentClock - wardStateTimes[targetWard];
        blockedTimes[targetWard].push_back(sample);
        wardStateTimes[targetWard] = currentClock;
        cout << "added blocked sample = " << sample << " to ward " << targetWard << endl;
    }
    
}


bool RelocSimulation::attemptDischarge(double & currentClock, int &serIdx, int &inService, vector<vector<int>> &wardOccupancy, vector<int> &capUse){
    
    currentClock = service_array[serIdx].serviceClock;
    
    //subtract patient from service
    inService--;
    
    //subtract patient from occupancy tracking vectors
    wardOccupancy[service_array[serIdx].wardTarget][service_array[serIdx].patientType]--;
    capUse[service_array[serIdx].wardTarget]--;
    
    updateServiceArray(serIdx,inService);
    
    cout << "attempted discharge success" << endl;
    
    return(true);
}

bool RelocSimulation::attemptAdmission(double &currentClock, int &arrIdx, vector<vector<int>> &wardOccupancy, vector<int> &capUse, int &inService){
    
    currentClock = arrival_array[arrIdx].arrivalClock;
    
    bool Ok = false;
    if (arrival_array[arrIdx].wardTarget==arrival_array[arrIdx].patientType &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }else if(capUse[arrival_array[arrIdx].patientType]==getWardCapacity(arrival_array[arrIdx].patientType) &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }
    
    if (Ok){
        //adjust occupancy
        wardOccupancy[arrival_array[arrIdx].wardTarget][arrival_array[arrIdx].patientType]++;
        capUse[arrival_array[arrIdx].wardTarget]++;
        
        //insert into service array
        service_array[inService] = Patient(arrival_array[arrIdx].arrivalClock,arrival_array[arrIdx].serviceTime,
                arrival_array[arrIdx].wardTarget,arrival_array[arrIdx].patientType);
        //calculate and insert clock at discharge
        service_array[inService].serviceClock = currentClock+arrival_array[arrIdx].serviceTime;
        
        //adjust number of patients currently in service
        inService++;
    }
    
    //move to next arrival
    arrIdx++;
    
    if (Ok){
        cout << "attempted admission success" << endl;
    }else{
        cout << "attempted admission failed" << endl;
    }
    
    return(Ok);
}
    
void RelocSimulation::updateServiceArray(int idx, int &inService){
    //removes a patient in service and adjusts the entire list
    //assumes the patient has been subtracted from inService
    
    if (idx<(inService-1)){
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
        int idx, mn = numeric_limits<double>::max();
        for (int i=0; i<inService; i++){
            if (service_array[i].serviceClock<mn){
               mn = service_array[i].serviceClock;
               idx = i;
            }
        }
    
        return(idx);
    }else{
        return(0);
    }
    
}

int RelocSimulation::minTimeSamples(vector<int> &openTimes, vector<int> &blockedTimes){
    
    int mn = numeric_limits<double>::max();
    
    for (int i=0; i<openTimes.size(); i++){
        if (mn<openTimes[i]){
            mn = openTimes[i];
        }
    }
    for (int i=0; i<blockedTimes.size(); i++){
        if (mn<blockedTimes[i]){
            mn = blockedTimes[i];
        }
    }
    
    return(mn);
}


double RelocSimulation::randomExponential(double rate){
    return(log(1-randomUniform(0,1))/(-rate));
}

double RelocSimulation::randomUniform(double from, double to){
    double r = (double)rand() / RAND_MAX;
    return(from + r * (to - from));
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