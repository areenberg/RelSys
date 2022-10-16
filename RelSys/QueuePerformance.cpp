/*
 * Copyright 2022 Anders Reenberg Andersen.
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
 * File:   QueuePerformance.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on October 15, 2022, 4:40 PM
 */

#include "QueuePerformance.h"

#include <vector>
#include <math.h>

using namespace std;

QueuePerformance::QueuePerformance(int nQ, QueueData * wda, vector<double> tPoints):
nQueues(nQ),
wd_array(wda),
timePoints(tPoints)        
{
    initialize();
}

QueuePerformance::QueuePerformance(const QueuePerformance& orig) {
}

QueuePerformance::~QueuePerformance() {
}

void QueuePerformance::initialize(){
 
    //current occupancy distributions
    occupancy.resize(nQueues,0);
    
    //tracking of occupancy distributions during simulation
    occupancyFreq.resize(nQueues);
    for (int widx; widx<nQueues; widx++){
        occupancyFreq[widx].resize(timePoints.size());
        for (int tidx=0; tidx<timePoints.size(); tidx++){
            occupancyFreq[widx][tidx].resize( ((wd_array+widx)->capacity + 1) ,0);
        }
    }
    
    lastClock=0;
    cycleLen = wd_array->timeDep.size(); //get length of one cycle
}



void QueuePerformance::arrival(double &newClock, int &widx, bool track){
    
    //update tracking
    if (track){
        trackOccupancy(newClock,widx);
    }
    
    //update current occupancy
    occupancy[widx]++;
    lastClock=newClock;
}

void QueuePerformance::discharge(double &newClock, int &widx, bool track){
    
    //update tracking
    if (track){
        trackOccupancy(newClock,widx);
    }
    
    //update current occupancy
    occupancy[widx]--;
    lastClock=newClock;
}


double QueuePerformance::lastTimePoint(){
    //return the time point corresponding to the
    //clock in cl
    return(lastClock-floor(lastClock/cycleLen)*cycleLen);
}

void QueuePerformance::trackOccupancy(double &newClock, int &widx){
        
        delta=newClock-lastClock;
        tp=lastTimePoint(); //initialize origin
        tpoffset=0.0;
        tpidx=0;
        
        while (delta>0){
            
            while ((timePoints[tpidx]+tpoffset)<=tp){
                tpidx++;
                if (tpidx==timePoints.size()){
                    tpidx=0;
                    tpoffset+=cycleLen;
                }
            }
            delta-=((timePoints[tpidx]+tpoffset)-tp);
            if (delta>0.0){
                //record occupancy at time point index
                occupancyFreq[widx][tpidx][occupancy[widx]]++;
                
                tp=timePoints[tpidx]; //move origin
                tpoffset=0.0;
            }
        }    
}

vector<vector<double>> QueuePerformance::getWardDenDist(int widx){
    
    double sm;
    vector<vector<double>> denDist(timePoints.size());
    for (int tidx=0; tidx<timePoints.size(); tidx++){
        denDist[tidx].resize( ((wd_array+widx)->capacity+1) ,0);
        
        sm=0.0;
        for (int j=0; j<denDist[tidx].size(); j++){
            sm+=occupancyFreq[widx][tidx][j];
        }
        for (int j=0; j<denDist[tidx].size(); j++){
            denDist[tidx][j] = (double)occupancyFreq[widx][tidx][j]/sm;
        }
    }
    
    return(denDist);
}
    





