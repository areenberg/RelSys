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
        vector<int> prefQ,vector<int> evalQ, string mdlt):
modelType(mdlt),
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
    
}

void Model::setDefaultParameters(){
    
    
    
    
}

void Model::setupCustomers(){
    
    
}

void Model::setupQueues(){
    
    
}

void Model::determineModelType(){
    
    
    
}


void Model::runHeuristic(int main_widx){
    
    
    
}


void Model::runSimulation(){
    
    
    
}


void Model::runModel(){
    
    if (modelType=="auto"){
        determineModelType();
    }
    
    
}



double Model::estimateRuntime(int stateSpaceSize){
    //estimates the runtime of the approximation
    //in milliseconds using the size of the state space
    
    //runtime (milliseconds) = a*stateSpaceSize^2 + b*stateSpaceSize
    double a=2.309390e-08;
    double b=6.279563e-02;
    
    return((a*pow((double)stateSpaceSize,2.0)+b*(double)stateSpaceSize));
}