
/* 
 * File:   HyberQueue.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on September 15, 2020, 11:15 AM
 */

#include "HyperQueue.h"
#include "RelocSimulation.h"
#include "PhaseFitter.h"

#include <vector>
#include <iostream>
#include <valarray>
#include <limits.h>

using namespace std;

HyperQueue::HyperQueue(int widx, int statesBlocked, int statesOpen, double aRate, double sRate, RelocSimulation *sm):
wardIndex(widx),
openRates(statesOpen,0),
openDist(statesOpen,0),
blockedRates(statesBlocked,0),
blockedDist(statesBlocked,0),        
//number of states
numberOfStates(statesBlocked + statesOpen),
//queue parameters
arrivalRate(aRate),
serviceRate(sRate),
//simulation        
sim_pointer(sm)
{
}

HyperQueue::HyperQueue(const HyperQueue& orig) {
}

HyperQueue::~HyperQueue() {
}


void HyperQueue::fitOpenPH(int seed){
    //fits the parameters for the PH
    //distribution accounting for the open rates


    if (!sim_pointer->openTimes[wardIndex].empty()){
    
        PhaseFitter ph_open;
        int EMiterations = 100;
    
        ph_open.setInputSample(sim_pointer->openTimes[wardIndex]);
    
        //set output PH distribution
        ph_open.setHyberExponential(openRates.size());
    
        //run calculations
        ph_open.run(EMiterations,seed);
    
        //get the result
        //cout << "Exit-rates:" << endl;
        for (int i=0; i<openRates.size(); i++){
            //cout << ph_open.exit_rate_vector[i] << endl;
            openRates[i] = ph_open.exit_rate_vector[i];
        }
        //cout << "Distribution, pi:" << endl;
        for (int i=0; i<openDist.size(); i++){
            //cout << ph_open.init_dist[i] << endl;
            openDist[i] = ph_open.init_dist[i];
        }
    
    }else{
        lowLoadOpenPH();
    }
    
}

void HyperQueue::fitBlockedPH(int seed){
    //fits the parameters for the PH
    //distribution accounting for the closed rates
    
    
    if (!sim_pointer->blockedTimes[wardIndex].empty()){
    
        PhaseFitter ph_blocked;
        int EMiterations = 100;
    
        ph_blocked.setInputSample(sim_pointer->blockedTimes[wardIndex]);
    
        //set output PH distribution
        ph_blocked.setHyberExponential(blockedRates.size());
    
        //run calculations
        ph_blocked.run(EMiterations,seed);
    
        //get the result
        for (int i=0; i<blockedRates.size(); i++){
            blockedRates[i] = ph_blocked.exit_rate_vector[i];
        }
        for (int i=0; i<blockedDist.size(); i++){
            blockedDist[i] = ph_blocked.init_dist[i];
        }
    
    }else{
        lowLoadBlockedPH();
    }
    
    
}

void HyperQueue::fitAll(int seed){
    //fits all PH parameters
    
    //cout << "WARD " << wardIndex << endl;
    cout << "Open Parameters (phases=" << openRates.size() << ")" << endl;
    fitOpenPH(seed);
    cout << "Blocked Parameters (phases=" << blockedRates.size() << ")" << endl; 
    fitBlockedPH(seed);
    
}

void HyperQueue::lowLoadOpenPH(){
    //uses a heuristic approach to estimating
    //the PH parameters - used when the load of
    //a ward is very low
    
    for (int i=0; i<openRates.size(); i++){
        openRates[i] = 1e-16;
    }
    for (int i=0; i<openDist.size(); i++){
        openDist[i] = 1.0/(double)openDist.size();
    }
    
}
        
void HyperQueue::lowLoadBlockedPH(){
    //uses a heuristic approach to estimating
    //the PH parameters - used when the load of
    //a ward is very low
    
    for (int i=0; i<blockedRates.size(); i++){
        blockedRates[i] = 1e16;
    }
    for (int i=0; i<blockedDist.size(); i++){
        blockedDist[i] = 1.0/(double)blockedDist.size();
    }
    
}        