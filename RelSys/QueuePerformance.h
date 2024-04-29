
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
#include <string>

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
    
    //get the results
    vector<vector<double>> getWardDenDist(int widx); //return the results as relative frequencies
    vector<vector<double>> getMeanOccupancy(int widx); //return the mean occupancy including 95% conf. intervals
    double expectedOccupancy(vector<double> dist);
    double stanDevOccupancy(double mean, vector<double> dist);
    int countSamples(int widx, int tidx);
    void saveResults(string fileName, int widx);
    
    
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







