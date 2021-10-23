/*
 * Copyright 2020 Anders Reenberg Andersen.
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
    
}

void HyperQueue::fitBlockedPH(int seed){
    //fits the parameters for the PH
    //distribution accounting for the closed rates
    
    PhaseFitter ph_blocked;
    int EMiterations = 100;
    
    ph_blocked.setInputSample(sim_pointer->blockedTimes[wardIndex]);
    
    //set output PH distribution
    ph_blocked.setHyberExponential(blockedRates.size());
    
    //run calculations
    ph_blocked.run(EMiterations,seed);
    
    //get the result
//    cout << "Rates:" << endl;
    for (int i=0; i<blockedRates.size(); i++){
//        cout << ph_blocked.exit_rate_vector[i] << endl;
        blockedRates[i] = ph_blocked.exit_rate_vector[i];
    }
//    cout << "Distribution:" << endl;
    for (int i=0; i<blockedDist.size(); i++){
//        cout << ph_blocked.init_dist[i] << endl;
        blockedDist[i] = ph_blocked.init_dist[i];
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
