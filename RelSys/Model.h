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
 * File:   Model.h
 * Author: Anders Reenberg Andersen
 *
 * Created on 9. februar 2023, 15.44
 */

//Straight-forward evaluation of the system comprising
//both the RelocEvaluation (approximation method) and
//the RelocSimulation (simulation method)

#ifndef MODEL_H
#define MODEL_H

#include "Model.h"
#include "CustomerData.h"
#include "QueueData.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "SystemParameters.h"
#include "QueuePerformance.h"

#include <vector>
#include <string>

class Model {
public:
    
    
    //METHODS
    
    //evaluate the model and store the results
    //will automatically select the model type if modelType=auto
    void runModel();
    
    void setSeed();
    void setBurnIn();
    void setMinimumSimulationTime();
    void setMinSamples();
    void setHyperStates(int openStates, int blockedStates);
    
    //VARIABLES
    vector<vector<double>> queueDenDist;
    vector<vector<int>> queueFreqDist;
    vector<double> blockingProbability,
    expectedOccupancy,expOccFraction;
    
    
    Model(vector<double> arrRates, vector<double> serTimes,
        vector<int> cap,vector<vector<double>> relProbMat,
        vector<int> prefQ,vector<int> evalQ, string mdlt="auto");
    
    Model(const Model& orig);
    virtual ~Model();
private:

    
    //METHODS
    void initialize();
    void setDefaultParameters();
    void determineModelType();
    void runHeuristic(int main_widx);
    void runSimulation();
    void setupCustomers();
    void setupQueues();
    double estimateRuntime(int stateSpaceSize);
    
    
    //VARIABLES
    RelocEvaluation * mdlHeu; //pointer for the approximation
    RelocSimulation * mdlSim; //pointer for the simulation
    QueueData * queues;
    CustomerData * customers; 
    
    string modelType; //indicate the model chosen (auto, simulation or approximation)
    vector<int> evalQueues; //indices of queues to evaluate from the QueueData array
    int nQueues, seed, minSamples, 
        openHyperStates, blockedHyperStates;
    double burnIn, minTime;
    
    //fundamental system parameters
    vector<double> arrivalRates,serviceTimes;
    vector<int> capacity, preferredQueues, evaluateQueues;
    vector<vector<double>> relocationProbabilityMatrix;
    
    
};

#endif /* MODEL_H */






