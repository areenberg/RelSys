
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
#include "QueuePerformance.h"

#include <iostream>
#include <vector>
#include <string>

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
    
    //time-dependent arrival rate as a fraction of the maximal arrival rate for the queue
    //note: only applicable to the simulation module
    //<<UNCOMMENT BELOW TO APPLY TIME-DEPENDENCY>>
//    vector<vector<double>> timeDep = {{0.5,0.8,0.9,0.01,0.5,1.0,1.0},
//                                      {0.6,0.7,1.0,0.9,0.7,0.5,0.4},
//                                      {0.1,0.1,0.2,0.3,0.6,1.0,0.8},
//                                      {0.5,1.0,1.0,1.0,0.5,0.3,0.1}};
    
    
    //now create the queue objects
    QueueData * wd_array = new QueueData[nQueues];
    for (int i=0; i<nQueues; i++){
        wd_array[i] = QueueData(i, //<<-- important: this input must always equal the index of the queue in the array
                                sysParam.queueArrivalRate(i),
                                sysParam.queueServiceRate(i),
                                capacity[i],
                                sysParam.queueRelProbability(i));
//        wd_array[i].addTimeDependency(timeDep[i]); //<<-- UNCOMMENT TO APPLY TIME-DEPENDENCY
//        wd_array[i].equalServiceRate(); //<<-- UNCOMMENT TO EQUALIZE SERVICE RATES
    }

    //--------------------------
    //CTMC APPROX EVALUATION
    //--------------------------

    //setup model object
    RelocEvaluation mdl(nQueues,wd_array);
    
    //first simulate open/blocked time-windows
    int seed = 123;
    int bin = -1; //set to -1 for auto
    int mt = 1; //set to 1 for auto
    mdl.runSimulation(seed,bin,mt,50);
    
    //choose a queue to evaluate
    int widx = 0; //queue index to be evaluated
    
    //choose number of states in PH distributions (optional)
    mdl.setOpenHyperStates(2);
    mdl.setBlockedHyperStates(1);
    
    //run the model
    mdl.runHeuristic(widx);
    
    //mdl.validateModel(widx);
    //mdl.evaluateModel();
    
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
    double burnIn = -1; //burn-in time (set to -1 for auto)
    double minTime = -1; //minimum simulation time (set to -1 for auto)
    vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
    sim_mdl.disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling
    
    //<<UNCOMMENT BELOW TO APPLY TIME-DEPENDENCY>>
//    vector<double> timePoints = {0.5,1.5,2.5,3.5,4.5,5.5,6.5}; //choose time points to track
//    QueuePerformance * qPer = new QueuePerformance(nQueues,wd_array,timePoints); //use this object to get the results
//    sim_mdl.enableTimeDependency(qPer); //enable time-dependency in simulation module
    
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

    //<<UNCOMMENT BELOW TO GET RESULTS FOR TIME-DEPENDENT SYSTEM>>
//    cout << "Time dep. occupancy:" << endl;
//    vector<vector<double>> dd = qPer->getWardDenDist(sim_widx);
//    for (int tidx=0; tidx<dd.size(); tidx++){
//        cout << timePoints[tidx] << ": ";
//        for (int i=0; i<dd[tidx].size(); i++){
//            cout << dd[tidx][i] << " ";
//        }
//        cout << endl;
//    }
//    //save results for all queues
//    for (int widx=0; widx<nQueues; widx++){
//        string fileName = "time_dep_occupancy";
//        fileName.append(to_string(widx));
//        fileName.append(".csv");
//        qPer->saveResults(fileName,widx);
//    }
    
    return 0;
}