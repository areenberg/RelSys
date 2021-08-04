
/* 
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on August 4, 2021, 02:11 PM
 */

#include "HeuristicQueue.h"
#include "HyperQueue.h"
#include "LinSolver.h"
#include "WardData.h"
#include "EntireSystem.h"
#include "RelocSimulation.h"
#include "PhaseFitter.h"
#include "RelocEvaluation.h"

#include <iostream>
#include <vector>

using namespace std;


int main(int argc, char** argv) {

    //--------------------------
    //QUEUE SETUP
    //--------------------------
    
    //total number of queues    
    int nWards = 4; 
    
    //arrival rates for each queue
    vector<double> arrivalRates = {1.0,0.5,2.5,2.0};
    
    //service rates for each queue
    vector<double> serviceRates = {0.04,0.02,0.1,0.1};
    
    //capacity of each queue
    vector<int> capacity = {30,40,10,15};
    
    //fraction of rejected customers that are moved to an alternative queue node
    vector<vector<double>> relProbs = {{0.0,0.1,0,8,0.1},
                                       {0.1,0.0,0.3,0.1}, //<<-- note: these do not have to sum to one
                                       {0.4,0.5,0.0,0.1},
                                       {0.2,0.5,0.3,0.0}};
    
    //now create the queue objects
    WardData * wd_array = new WardData[nWards];
    wd_array[0] = WardData(0,arrivalRates[0],serviceRates[0],capacity[0],relProbs[0]);
    wd_array[1] = WardData(1,arrivalRates[1],serviceRates[1],capacity[1],relProbs[1]);
    wd_array[2] = WardData(2,arrivalRates[2],serviceRates[2],capacity[2],relProbs[2]);
    wd_array[3] = WardData(3,arrivalRates[3],serviceRates[3],capacity[3],relProbs[3]);
    
    //--------------------------
    //HEURISTIC EVALUATION
    //--------------------------
    
    //setup model object
    RelocEvaluation mdl(nWards,wd_array);
    
    //first simulate open/blocked time-windows
    int seed = 123;
    mdl.runSimulation(seed,365,365,50);
    
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
    RelocSimulation sim_mdl(nWards,wd_array);
        
    //setup and run simulation
    sim_mdl.setSeed(123); //set the seed
    double burnIn = 365; //burn-in time
    double minTime = 50000; //minimum simulation time
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

