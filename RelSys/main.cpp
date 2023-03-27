
/*
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 */

#include "Model.h"

#include <iostream>
#include <vector>

using namespace std;


//--------------------------
//STRUCTURES
//--------------------------

//input data structure (including example values)
struct {
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};
    vector<double> serviceTimes {10,5,10,8};
    vector<int> capacity = {15,20,10,30};
    vector<vector<double>> relocationProbabilities = {{0.0,0.4,0.1,0.5},
                                                      {0.3,0.0,0.5,0.0},
                                                      {0.0,0.5,0.0,0.5},
                                                      {0.2,0.3,0.5,0.0}};
    vector<int> preferredQueue = {0,1,2,3};
} data;

//model settings structure
struct {
    bool verbose=false;
    bool evalAll=true;
    bool equalizeService=false;
    vector<int> evaluatedQueue;
    string modelType="simulation";

    int seed=-1;
    double burnIn=-1;
    double minimumSimulationTime=-1;
    int minSamples=-1;
    vector<int> hyperStates = {-1,-1};
    double accTol=5e-3;
    string accSampleType="preferred";
} settings;

//results
struct {
    bool evaluated=false;
    vector<double> shortageProbability,shortageProbabilityPref;
    vector<double> availProbability,availProbabilityPref;
    vector<vector<double>> queueDenDist,queueDenDistPref;
    vector<vector<int>> queueFreqDist,queueFreqDistPref;
    vector<double> expectedOccupancy;
    vector<double> expOccFraction;
} results;


//--------------------------
//METHODS
//--------------------------

// READ ARGUMENTS



//OTHER METHODS

void setVerbose(bool set){
    settings.verbose=set;
}

void evaluateAllQueues(){
    settings.evaluatedQueue.resize(data.capacity.size());
    for (int i=0; i<data.capacity.size(); i++){
        settings.evaluatedQueue[i] = i;
    }
}

void queuesToEvaluate(vector<int> qEvalIdx){
    //set the indices of queues to
    //evaluate
    settings.evalAll=false;
    settings.evaluatedQueue = qEvalIdx;
}

void setType(string mdltype){
    //set the type of model to use
    //in the evaluation
    settings.modelType = mdltype;
}

void equalizeService(bool equalize){
    //specifies if the service rates
    //in the model should be equalized
    //and loads correspondingly balanced
    settings.equalizeService = equalize;
}

void setAccuracySampleType(string stype){
    settings.accSampleType = stype;    
}

void setSimulationTolerance(double tol){
    settings.accTol = tol;
}

void setSeed(int sd){
    settings.seed = sd;
}

void setBurnIn(double bn){
    if (bn>0){
        settings.burnIn=bn;
    }else{
        cout << "Burn-in time must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setMinimumSimulationTime(double mnTime){
    if (mnTime>0){
        settings.minimumSimulationTime=mnTime;
    }else{
        cout << "Simulation time must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setMinSamples(int mnSamples){
    if (mnSamples>0){
        settings.minSamples=mnSamples;
    }else{
        cout << "Min. number of open/shortage samples must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setHyperStates(int openStates, int blockedStates){
    if (openStates>0 && blockedStates>0){
        settings.hyperStates[0] = openStates;
        settings.hyperStates[1] = blockedStates;
    }else{
        cout << "The number of open and blocked states must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void checkParameters(){
    //checks if parameters are feasible

    if (data.arrivalRates.size() != data.serviceTimes.size()){
        cout << "The number of arrival rates must equal the number of service times. Aborting program." << endl;
        exit(1);
    }
    if (data.arrivalRates.size() != data.relocationProbabilities.size()){
        cout << "The number of arrival rates, and service times, must equal the number of rows in the relocation matrix. Aborting program." << endl;
        exit(1);
    }
    if (data.arrivalRates.size() != data.preferredQueue.size()){
        cout << "The number of arrival rates, and service times, must equal the length of the vector specifying the preferred queues. Aborting program." << endl;
        exit(1);
    }
    if (settings.evaluatedQueue.size()>data.capacity.size() || settings.evaluatedQueue.empty()){
        cout << "The length of the vector of evaluated queues must be between 1 and n_queues. Aborting program." << endl;
        exit(1);
    }        
    for (int i=0; i<data.relocationProbabilities.size(); i++){
        if (data.capacity.size() != data.relocationProbabilities[i].size()){
            cout << "The length of the capacity vector must equal the number of columns in the relocation matrix. Aborting program." << endl;
            exit(1);
        }
    }        
    for (int i=0; i<data.arrivalRates.size(); i++){
        if (data.arrivalRates[i]<=0.0 || data.serviceTimes[i]<=0.0 ){
            cout << "Arrival rates and service times must be larger than 0.0. Aborting program." << endl;
            exit(1);
        }
        for (int j=0; j<data.relocationProbabilities[i][j]; j++){
            if (data.relocationProbabilities[i][j]<0.0 || data.relocationProbabilities[i][j]>1.0){
                cout << "Values in the relocation matrix must be between 0.0 and 1.0. Aborting program." << endl;
                exit(1);
            }
        }
    }
    for (int i=0; i<data.preferredQueue.size(); i++){
        if (data.preferredQueue[i]<0 || data.preferredQueue[i]>(data.capacity.size()-1)){
            cout << "Indices in the vector of preferred queues must be between 0 and n_queues-1. Aborting program." << endl;
            exit(1);
        }
    }
    for (int i=0; i<data.capacity.size(); i++){
        if (data.capacity[i]<1){
            cout << "The capacity of each queue must be equal to or larger than 1. Aborting program." << endl;
            exit(1);
        }
    }          
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        if (settings.evaluatedQueue[i]<0 || settings.evaluatedQueue[i]>(data.capacity.size()-1)){
            cout << "Indices in the vector of evaluated queues must lie in the interval between 0 and n_queues-1. Aborting program." << endl;
            exit(1);
        }    
    }
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        for (int j=0; j<settings.evaluatedQueue.size(); j++){
            if (i!=j && settings.evaluatedQueue[i]==settings.evaluatedQueue[j]){
                cout << "All indices in the vector of evaluated queues must be unique. The same queue cannot be evaluated twice. Aborting program." << endl;
                exit(1);
            }
        }    
    }
    double sm;
    for (int i=0; i<data.relocationProbabilities.size(); i++){
        sm=0;
        for (int j=0; j<data.relocationProbabilities[i].size(); j++){
            sm+=data.relocationProbabilities[i][j];
        }
        if (sm>1.0){
            string out = "The sum of the relocation probabilities in row "+to_string(i)+" is equal to "+to_string(sm)+". The sum must be equal to or smaller than 1.0. Aborting program.";
            cout << out << endl;
            exit(1);
        }
    }
    if ((settings.minimumSimulationTime!=-1 && settings.burnIn==-1) || (settings.minimumSimulationTime==-1 && settings.burnIn!=-1)){
        cout << "It is not possible to set the overall simulation time without setting the burn-in time and vice versa. Aborting program." << endl;
        exit(1);
    }
    if (settings.minimumSimulationTime!=-1 && settings.burnIn!=-1 && settings.minimumSimulationTime<=settings.burnIn){
        cout << "The simulation time has to be longer than the burn-in time. Aborting program." << endl;
        exit(1);
    }
    if (settings.accTol<=0.0 || settings.accTol>=1.0){
        cout << "The tolerance level for the accuracy estimation procedure must be between 0 and 1 (e.g. 1e-2). Aborting program." << endl;
        exit(1);
    }        

}

void runCalculations(){

    //--------------------------
    //MODEL DATA
    //--------------------------

    //indices of queues to be evaluated by the model
    if (settings.evalAll){
        evaluateAllQueues();        
    }

    //check parameters are feasible
    checkParameters();    
    
    //--------------------------
    //MODEL EVALUATION
    //--------------------------

    //control verbose from model (deactivate cout)
    if (!settings.verbose){
        cout.setstate(std::ios_base::failbit);
    }

    //create the model object
    Model mdl(data.arrivalRates,data.serviceTimes,data.capacity,
            data.relocationProbabilities,data.preferredQueue,
            settings.evaluatedQueue,settings.modelType,
            settings.equalizeService);
    
    //additional settings
    if (settings.seed!=-1){
        mdl.setSeed(settings.seed);
    }
    if (settings.burnIn!=-1){
        mdl.setBurnIn(settings.burnIn);
    }
    if (settings.minimumSimulationTime!=-1){
        mdl.setMinimumSimulationTime(settings.minimumSimulationTime);
    }
    if (settings.minSamples!=-1){
        mdl.setMinSamples(settings.minSamples);
    }
    if (settings.hyperStates[0]!=-1 && settings.hyperStates[1]!=-1){
        mdl.setHyperStates(settings.hyperStates[0],settings.hyperStates[1]);
    }
    mdl.setSimTolerance(settings.accTol);
    mdl.setAccuracySampleType(settings.accSampleType);
    
    //now evaluate the model
    mdl.runModel();

    cout.clear(); //reactivate cout
    
    //--------------------------
    // GET THE RESULTS
    //--------------------------

    //queue density distribution    
    results.queueDenDist.resize(data.capacity.size());
    results.queueDenDistPref.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueDenDist[settings.evaluatedQueue[i]].resize(mdl.queueDenDist[settings.evaluatedQueue[i]].size());
        results.queueDenDistPref[settings.evaluatedQueue[i]].resize(mdl.queueDenDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueDenDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueDenDist[settings.evaluatedQueue[i]][j] = mdl.queueDenDist[settings.evaluatedQueue[i]][j];
            results.queueDenDistPref[settings.evaluatedQueue[i]][j] = mdl.queueDenDistPref[settings.evaluatedQueue[i]][j];
        }
    }

    //queue frequency distribution    
    results.queueFreqDist.resize(data.capacity.size());
    results.queueFreqDistPref.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueFreqDist[settings.evaluatedQueue[i]].resize(mdl.queueFreqDist[settings.evaluatedQueue[i]].size());
        results.queueFreqDistPref[settings.evaluatedQueue[i]].resize(mdl.queueFreqDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueFreqDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueFreqDist[settings.evaluatedQueue[i]][j] = mdl.queueFreqDist[settings.evaluatedQueue[i]][j];
            results.queueFreqDistPref[settings.evaluatedQueue[i]][j] = mdl.queueFreqDistPref[settings.evaluatedQueue[i]][j];
        }
    }
    
    //shortage probabilities
    results.shortageProbability.resize(data.capacity.size(),-1);
    results.shortageProbabilityPref.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.shortageProbability[settings.evaluatedQueue[i]] = mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.shortageProbabilityPref[settings.evaluatedQueue[i]] = mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
    }

    //availability probabilities
    results.availProbability.resize(data.capacity.size(),-1);
    results.availProbabilityPref.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.availProbability[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.availProbabilityPref[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
    }

    //expected occupancy
    results.expectedOccupancy.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expectedOccupancy[settings.evaluatedQueue[i]] = mdl.expectedOccupancy[settings.evaluatedQueue[i]];
    }

    //expected fraction of customers occupied
    results.expOccFraction.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expOccFraction[settings.evaluatedQueue[i]] = mdl.expOccFraction[settings.evaluatedQueue[i]];
    }
    
    results.evaluated = true;

}

int main(int argc, char** argv) {

    //--------------------------
    // READ ARGUMENTS
    //--------------------------
    
    //CODE HERE
    
    //--------------------------
    // EVALUATE MODEL
    //--------------------------
    
    runCalculations();

    //--------------------------
    // DISPLAY OR SAVE RESULTS
    //--------------------------
    
    //CODE HERE
    

    return 0;
}