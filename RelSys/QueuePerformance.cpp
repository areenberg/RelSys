
/* 
 * File:   QueuePerformance.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on October 15, 2022, 4:40 PM
 */

#include "QueuePerformance.h"

#include <vector>
#include <math.h>
#include <fstream>

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

vector<vector<double>> QueuePerformance::getMeanOccupancy(int widx){
    
    vector<vector<double>> dd = getWardDenDist(widx);
    double stanDev; int n;
    
    vector<vector<double>> res(timePoints.size());
    for (int tidx=0; tidx<timePoints.size(); tidx++){
        res[tidx].resize(3,0); //lower conf., mean, upper conf.
        res[tidx][1] = expectedOccupancy(dd[tidx]);
        stanDev=stanDevOccupancy(res[tidx][1],dd[tidx]);
        n=countSamples(widx,tidx);
        res[tidx][0] = res[tidx][1] - (1.959964 * stanDev/sqrt((double)n));
        res[tidx][2] = res[tidx][1] + (1.959964 * stanDev/sqrt((double)n));
    }
    
    return(res);
}

int QueuePerformance::countSamples(int widx, int tidx){
    
    int sm=0.0;
    for (int j=0; j<occupancyFreq[widx][tidx].size(); j++){
        sm+=occupancyFreq[widx][tidx][j];
    }
    
    return(sm);
}


double QueuePerformance::expectedOccupancy(vector<double> dist){
    
    double y=0.0;
    for (int i=0; i<dist.size(); i++){
        y += i*dist[i];
    }
    
    return(y);
}
    
double QueuePerformance::stanDevOccupancy(double mean, vector<double> dist){
    
    double y=0.0;
    for (int i=0; i<dist.size(); i++){
        y += pow((i-mean),2)*dist[i];
    }
    y = sqrt(y);
    
    return(y);
}

void QueuePerformance::saveResults(string fileName, int widx){
    
    ofstream resFile(fileName);
    
    vector<vector<double>> dd = getWardDenDist(widx);
    
    //create the header
    resFile << "time_point,";
    for (int i=0; i<(wd_array+widx)->capacity; i++){
        resFile << i << ",";
    }
    resFile << (wd_array+widx)->capacity << endl;
    
    for (int tidx=0; tidx<dd.size(); tidx++){
        resFile << timePoints[tidx] << ",";
        for (int i=0; i<(dd[tidx].size()-1); i++){
            resFile << dd[tidx][i] << ",";
        }
        resFile << dd[tidx][(dd[tidx].size()-1)] << endl;
    }    
    resFile.close();
    
}




