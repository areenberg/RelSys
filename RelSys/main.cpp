
/*
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on August 4, 2021, 02:11 PM
 */

#include "CustomerData.h"
#include "QueueData.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "SystemParameters.h"

#include <iostream>
#include <vector>

using namespace std;


int main(int argc, char** argv) {

    //--------------------------
    //CUSTOMER SETUP
    //--------------------------

    //total number of customer types
    int nCustomerTypes = 4;

    //arrival rates for each customer type
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};

    //mean service time for each customer type
    vector<double> serviceTimes = {10,5,10,8};

    //fraction of rejected customers that are moved to an alternative queue node
    //this is an nCustomerTypes x nQueues matrix
    vector<vector<double>> relProbs = {{0.0,0.4,0.1,0.5},
                                       {0.3,0.0,0.5,0.0}, //<<-- note: these do not have to sum to one
                                       {0.0,0.5,0.0,0.5},
                                       {0.2,0.3,0.5,0.0}};
    
    //queue indices preferred by each customer type
    vector<int> preferredQueue = {0,1,2,3};

    //create the customer type objects
    CustomerData * custs_array = new CustomerData[nCustomerTypes];
    for (int i=0; i<nCustomerTypes; i++){
        custs_array[i] = CustomerData(preferredQueue[i],
                                      arrivalRates[i],
                                      serviceTimes[i],
                                      relProbs[i]);
    }

    //--------------------------
    //QUEUE SETUP
    //--------------------------

    //total number of queues
    int nQueues = 4;

    //capacity of each queue
    vector<int> capacity = {15,20,10,30};

    //calculate system input parameters from customer types to queues
    SystemParameters sysParam(nQueues,nCustomerTypes,custs_array);

    //now create the queue objects
    QueueData * wd_array = new QueueData[nQueues];
    for (int i=0; i<nQueues; i++){
        wd_array[i] = QueueData(i, //<<-- important: this input must always equal the index of the queue in the array
                                sysParam.queueArrivalRate(i),
                                sysParam.queueServiceRate(i),
                                capacity[i],
                                sysParam.queueRelProbability(i));
    }

    //--------------------------
    //CTMC APPROX EVALUATION
    //--------------------------

    //setup model object
    RelocEvaluation mdl(nQueues,wd_array);

    //first simulate open/blocked time-windows
    int seed = 123;
    int bin = -1;
    int mt = 1;
    mdl.runSimulation(seed,bin,mt,50);
    
    //choose a queue to evaluate
    int widx = 0; //queue index to be evaluated
    mdl.runHeuristic(widx);

    //get the result
    cout << "--- RESULTS ---" << endl;

    cout << "Marginal probability distribution:" << endl;
    for (int i=0; i<mdl.marginalDist.size(); i++){
        cout << mdl.marginalDist[i] << endl;
    }

    cout << "Probability of rejection:" << endl;
    cout << mdl.blockingProbability << endl;

    cout << "Expected server occupancy:" << endl;
    cout << mdl.expectedOccupancy << endl;

    cout << "Expected fraction of servers occupied:" << endl;
    cout << mdl.expOccFraction << endl;

    //--------------------------
    //SIMULATION EVALUATION
    //--------------------------

    //setup simulation model object
    RelocSimulation sim_mdl(nQueues,wd_array);

    //setup and run simulation
    sim_mdl.setSeed(123); //set the seed
    double burnIn = -1; //burn-in time
    double minTime = -1; //minimum simulation time
    vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
    sim_mdl.disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling

    sim_mdl.simulate(burnIn,minTime,maxWardSamples); //now run the simulation

    //In contrary to the heuristic evaluation, the entire system is evaluated at the
    //same time. However, for the sake of this demonstration we choose only to print
    //the results from a single queue.
    int sim_widx = 0;

    //get the result
    cout << "--- RESULTS ---" << endl;

    cout << "Marginal frequency distribution:" << endl;
    for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
        cout << sim_mdl.wardFreqDist[sim_widx][i] << endl;
    }
    cout << "Marginal density distribution:" << endl;
    for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
        cout << sim_mdl.wardDenDist[sim_widx][i] << endl; //corresponds to marginalDist in the heuristic evaluation
    }

    cout << "Probability of rejection:" << endl;
    cout << sim_mdl.blockingProbability[sim_widx] << endl;

    cout << "Expected server occupancy:" << endl;
    cout << sim_mdl.expectedOccupancy[sim_widx] << endl;

    cout << "Expected fraction of servers occupied:" << endl;
    cout << sim_mdl.expOccFraction[sim_widx] << endl;

    return 0;
}
