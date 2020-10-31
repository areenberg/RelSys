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

#include <vector>
#include <iostream>

using namespace std;

HyperQueue::HyperQueue(int statesBlocked, int statesOpen, double aRate, double sRate):
openRates(statesOpen,0),
openDist(statesOpen,0),
blockedRates(statesBlocked,0),
blockedDist(statesBlocked,0),        
//number of states
numberOfStates(statesBlocked + statesOpen),
//queue parameters
arrivalRate(aRate),
serviceRate(sRate)
{    
}

HyperQueue::HyperQueue(const HyperQueue& orig) {
}

HyperQueue::~HyperQueue() {
}


void HyperQueue::fitOpenPH(int q){
    //fits the parameters for the PH
    //distribution accounting for the open rates
    
    //demo of parameter fitting
    if (q==0){
        openRates[0] = 0.809498; openRates[1] = 4.926081;
        openDist[0] = 0.348426; openDist[1] = 0.651574;
    }else if (q==1){
        openRates[0] = 0.042962; openRates[1] = 3.554449;
        openDist[0] = 0.426178; openDist[1] = 0.573822;
    }else if (q==2){
        openRates[0] = 0.246928; openRates[1] = 3.320117;
        openDist[0] = 0.383202; openDist[1] = 0.616798;
    }else if (q==3){
        openRates[0] = 0.083738; openRates[1] = 4.282059;
        openDist[0] = 0.386958; openDist[1] = 0.613042;
    }else{
        cout << "Wrong hyperqueue index." << endl;
    }
    
}

void HyperQueue::fitBlockedPH(int q){
    //fits the parameters for the PH
    //distribution accounting for the closed rates
    
    //demo of parameter fitting
    if (q==0){
        blockedRates[0] = 2.438961;
        blockedDist[0] = 1.0;
    }else if(q==1){
        blockedRates[0] = 6.776864;
        blockedDist[0] = 1.0;
    }else if(q==2){
        blockedRates[0] = 3.551759;
        blockedDist[0] = 1.0;        
    }else if(q==3){
        blockedRates[0] = 6.892878;
        blockedDist[0] = 1.0;    
    }else{
        cout << "Wrong hyperqueue index." << endl;
    }
    
    
}

void HyperQueue::fitAll(int q){
    //fits all PH parameters
    
    fitOpenPH(q);
    fitBlockedPH(q);
    
}
