
/*
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 */

#include "Model.h"

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

    //--------------------------
    //MODEL DATA
    //--------------------------
    //arrival rates for each customer type
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};
    
    //mean service time for each customer type
    vector<double> serviceTimes = {10,5,10,8};

    //capacity of each queue
    vector<int> capacity = {15,20,10,30};
    
    //fraction of rejected customers that are moved to an alternative queue node
    //this is an nCustomerTypes x nQueues matrix
    vector<vector<double>> relocationProbabilities = {{0.0,0.4,0.1,0.5},
                                                      {0.3,0.0,0.5,0.0}, //<<-- note: these do not have to sum to one
                                                      {0.0,0.5,0.0,0.5},
                                                      {0.2,0.3,0.5,0.0}};
    
    //queue indices preferred by each customer type
    vector<int> preferredQueue = {0,1,2,3};

    //indices of the queues to be evaluated    
    vector<int> evaluatedQueue = {0,1,2,3};

    //--------------------------
    //MODEL EVALUATION
    //--------------------------
    
    //create the model object
    Model mdl(arrivalRates,serviceTimes,capacity,
            relocationProbabilities,preferredQueue,
            evaluatedQueue);
    
    //now evaluate the model
    mdl.runModel();

    //--------------------------
    // GET THE RESULTS
    //--------------------------
    
    cout << endl << "#### RESULTS #####" << endl;

    cout << "Marginal density distribution:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << "----" << "Queue " << evaluatedQueue[i] << "----" << endl;
        for (int j=0; j<mdl.queueDenDistPref[evaluatedQueue[i]].size(); j++){
            cout << mdl.queueDenDistPref[evaluatedQueue[i]][j] << endl;
        }
    }
    
    cout << endl << "Probability of rejection:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.blockingProbabilityPref[evaluatedQueue[i]] << endl;
    }
    
    cout << endl << "Expected server occupancy:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.expectedOccupancy[evaluatedQueue[i]] << endl;
    }

    cout << endl << "Expected fraction of servers occupied:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.expOccFraction[evaluatedQueue[i]] << endl;
    }

    return 0;
}