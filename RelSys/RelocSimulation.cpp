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

void RelocSimulation::generateArrivalList(int length, double currentClock){
    
    double st;
    patient_array = new Patient[length];
    
    for (int i=0; i<length; i++){
        currentClock += randomExponential(1.0);
        st = randomExponential(0.5);
        patient_array[i] = Patient(currentClock,st,1,1);
    }
    
}

void RelocSimulation::simulate(double minTime, double minSamples){
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
    double clock;
    int patientArraySize = 10;
    generateArrivalList(patientArraySize,clock);
    
    for (int i=0; i<patientArraySize; i++){
        cout << patient_array[i].arrivalTime << " " << patient_array[i].serviceTime << endl;
    }
 
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