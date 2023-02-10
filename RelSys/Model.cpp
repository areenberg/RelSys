/*
 * Copyright 2023 Anders Reenberg Andersen.
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
 * File:   Model.cpp
 * Author: $Anders Reenberg Andersen
 * 
 * Created on 9. februar 2023, 15.44
 */

//Straight-forward evaluation of the system comprising
//both the RelocEvaluation (approximation method) and
//the RelocSimulation (simulation method)

#include "Model.h"
#include "CustomerData.h"
#include "QueueData.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "SystemParameters.h"
#include "QueuePerformance.h"

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>


using namespace std;

Model::Model(vector<double> arrRates, vector<double> serTimes,
        vector<int> cap,vector<vector<double>> relProbMat,
        vector<int> prefQ,vector<int> evalQ, string mdlt,
        bool eqze):
modelType(mdlt),
equalize(eqze),
arrivalRates(arrRates),
serviceTimes(serTimes),
capacity(cap), 
preferredQueues(prefQ),
evaluateQueues(evalQ),
relocationProbabilityMatrix(relProbMat)            
{
    initialize();
}

Model::Model(const Model& orig) {
}

Model::~Model() {
}



void Model::initialize(){
    
    setDefaultParameters();
    setupCustomers();
    setupQueues();
    prepareOutput();
    
}

void Model::setDefaultParameters(){
    
    nQueues = relocationProbabilityMatrix[0].size();
    nCustomers = arrivalRates.size();
    
    seed=123;
    minSamples=50;
    burnIn=-1;
    minTime=-1;
    openHyperStates=2; 
    blockedHyperStates=1;
    
}

void Model::setSeed(int sd){
    seed=sd;
}

void Model::setBurnIn(double bn){
    burnIn=bn;
}

void Model::setMinimumSimulationTime(double mnTime){
    minTime=mnTime;
}

void Model::setMinSamples(int mnSamples){
    minSamples=mnSamples;
}

void Model::setHyperStates(int openStates, int blockedStates){
    openHyperStates=openStates; 
    blockedHyperStates=blockedStates;
}

void Model::setupCustomers(){
    
    customers = new CustomerData[nCustomers];
    for (int i=0; i<nCustomers; i++){
        customers[i] = CustomerData(preferredQueues[i],
                                    arrivalRates[i],
                                    serviceTimes[i],
                                    relocationProbabilityMatrix[i]);
    }
    
}

void Model::setupQueues(){
    
    //consolidate customers into input flows of each queue
    sysParam = new SystemParameters(nQueues,nCustomers,customers);
    
    //create the queue objects
    queues = new QueueData[nQueues];
    for (int i=0; i<nQueues; i++){
        queues[i] = QueueData(i, //<<-- important: this input must always equal the index of the queue in the array
                                sysParam->queueArrivalRate(i),
                                sysParam->queueServiceRate(i),
                                capacity[i],
                                sysParam->queueRelProbability(i));
        if (equalize){
            queues[i].equalServiceRate(); //equalize the service rates
        }
    }
    
}

void Model::prepareOutput(){
    //prepare the vector containers
    //storing the results
    
    //prepares the vector for all queues
    //even though only some of them might
    //be evaluated
    
    queueDenDist.resize(nQueues);
    queueFreqDist.resize(nQueues);
    for (int i=0; i<nQueues; i++){
        queueDenDist[i].resize((capacity[i]+1),-1);
        queueFreqDist[i].resize((capacity[i]+1),-1);
    }
    blockingProbability.resize(nQueues,-1);
    expectedOccupancy.resize(nQueues,-1);
    expOccFraction.resize(nQueues,-1);
    
}

void Model::determineModelType(){
    
    
    
}


void Model::runHeuristic(int main_widx){
    //CTMC APPROX EVALUATION
    
    //setup model object
    mdlHeu = new RelocEvaluation(nQueues,queues);
    
    //first simulate open/blocked time-windows
    int mt;
    if (minTime<=0.0){
        mt=1;
    }else{
        mt=(int)minTime;
    }
    mdlHeu->runSimulation(seed,burnIn,mt,minSamples);
    
    //choose number of states in PH distributions (optional)
    mdlHeu->setOpenHyperStates(openHyperStates);
    mdlHeu->setBlockedHyperStates(blockedHyperStates);
    
    //run (validate and evaluate) the model
    mdlHeu->runHeuristic(main_widx);
    
    //get the result
    for (int i=0; i<mdlHeu->marginalDist.size(); i++){
        queueDenDist[main_widx][i] = mdlHeu->marginalDist[i];
    }
    blockingProbability[main_widx] = mdlHeu->blockingProbability;
    expectedOccupancy[main_widx] = mdlHeu->expectedOccupancy;
    expOccFraction[main_widx] = mdlHeu->expOccFraction;
    
    delete mdlHeu;
}


void Model::runSimulation(){
    //SIMULATION EVALUATION
    
    //setup simulation model object
    mdlSim = new RelocSimulation(nQueues,queues);
    
    //setup and run simulation
    mdlSim->setSeed(seed); //set the seed
    vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
    mdlSim->disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling
    mdlSim->simulate(burnIn,minTime,maxWardSamples); //now run the simulation


    //get the result
    int widx;
    for (int idx=0; idx<evaluateQueues.size(); idx++){
        widx=evaluateQueues[idx];
        
        for (int i=0; i<mdlSim->wardFreqDist[widx].size(); i++){
            queueFreqDist[widx][i] = mdlSim->wardFreqDist[widx][i];
        }
        for (int i=0; i<mdlSim->wardDenDist[widx].size(); i++){
            queueDenDist[widx][i] = mdlSim->wardDenDist[widx][i];
        }
        
        blockingProbability[widx] = mdlSim->blockingProbability[widx];
        expectedOccupancy[widx] = mdlSim->expectedOccupancy[widx];
        expOccFraction[widx] = mdlSim->expOccFraction[widx];
    }
    
    delete mdlSim;
}


void Model::runModel(){
    
    if (modelType=="auto"){
        determineModelType();
    }
    
    
}

int Model::evaluateStateSpaceSize(int main_widx){
    
    
    return(0);
}


double Model::estimateRuntime(int stateSpaceSize){
    //estimates the runtime of the approximation
    //in milliseconds using the size of the state space
    
    //runtime (milliseconds) = a*stateSpaceSize^2 + b*stateSpaceSize
    double a=2.309390e-08;
    double b=6.279563e-02;
    
    return((a*pow((double)stateSpaceSize,2.0)+b*(double)stateSpaceSize));
}