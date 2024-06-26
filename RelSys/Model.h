
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
    
    void setSeed(int sd);
    void setBurnIn(double bn);
    void setSimTolerance(double at);
    void setAccuracySampleType(string stype);
    void setMinimumSimulationTime(double mnTime);
    void setMinSamples(int mnSamples);
    void setHyperStates(int openStates, int blockedStates);
    
    //VARIABLES
    vector<vector<double>> queueDenDist,queueDenDistPref;
    vector<vector<int>> queueFreqDist,queueFreqDistPref;
    vector<double> blockingProbability,
    blockingProbabilityPref,expectedOccupancy,
    expOccFraction;
    
    
    Model(vector<double> arrRates, vector<double> serTimes,
        vector<int> cap,vector<vector<double>> relProbMat,
        vector<int> prefQ,vector<int> evalQ, string mdlt="simulation",
        bool eqze=false);
    
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
    void prepareOutput();
    double estimateRuntime(int stateSpaceSize);
    
    
    //VARIABLES
    RelocEvaluation * mdlHeu; //pointer for the approximation
    RelocSimulation * mdlSim; //pointer for the simulation
    QueueData * queues;
    CustomerData * customers;
    SystemParameters * sysParam;
    
    string modelType; //indicate the model chosen (auto, simulation or approximation)
    vector<int> evalQueues; //indices of queues to evaluate from the QueueData array
    int nQueues, nCustomers, seed, minSamples, 
        openHyperStates, blockedHyperStates, statSize;
    double accTol, burnIn, minTime;
    bool accPref,equalize; //indicates if service rates are equalized (default=true)
    
    //fundamental system parameters
    vector<double> arrivalRates,serviceTimes;
    vector<int> capacity, preferredQueues, evaluateQueues;
    vector<vector<double>> relocationProbabilityMatrix;
    
    
    
};

#endif /* MODEL_H */






