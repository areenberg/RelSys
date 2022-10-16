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
 * File:   QueuePerformance.h
 * Author: Anders Reenberg Andersen
 *
 * Created on October 15, 2022, 4:40 PM
 */

#ifndef QUEUEPERFORMANCE_H
#define QUEUEPERFORMANCE_H

#include "QueueData.h"

#include <vector>

using namespace std;

class QueuePerformance {
public:
    
    //occupancy distributions are sampled at these
    //time segments in the cycle
    vector<double> timePoints;
    
    //current occupancy distributions
    vector<int> occupancy; //ward index
    
    //sampled occupancy distributions
    vector<vector<vector<int>>> occupancyFreq; //ward index x time point index x occupancy 
    
    double lastClock;
    
    void arrival(double &newClock, int &widx, bool track=true);
    void discharge(double &newClock, int &widx, bool track=true);
    
    //return the results as relative frequencies
    vector<vector<double>> getWardDenDist(int widx);
    
    //dummy constructor (not included in cpp-file) 
    QueuePerformance() {};
    QueuePerformance(int nQ, QueueData * wda, vector<double> tPoints);
    QueuePerformance(const QueuePerformance& orig);
    virtual ~QueuePerformance();
private:
    
    int nQueues,cycleLen,tpidx;
    QueueData * wd_array;
    double delta,tp,tpoffset;
    
    void initialize();
    double lastTimePoint();
    void trackOccupancy(double &newClock, int &widx);
    
};

#endif /* QUEUEPERFORMANCE_H */







