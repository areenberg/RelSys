
/* 
 * File:   QueueData.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on November 1, 2020, 7:14 PM
 */

#include "QueueData.h"
#include "Combinatorial.h"

QueueData::QueueData(int wdnr, double arrRate, double serRate, int cap, vector<double> relProbs):
wardnr(wdnr),
arrivalRate(arrRate),
serviceRate(serRate),        
capacity(cap),
relocationProbabilities(relProbs)
{   
    //a value of -1 in the first index indicates that time-dependency is disabled
    timeDep.resize(1,-1.0);
}

QueueData::QueueData(const QueueData& orig) {
}

QueueData::~QueueData() {
}


void QueueData::equalServiceRate(){
    
    //use this method to set the service rate to a value of 1
    //and scale the arrival rate keeping the load on the system
    arrivalRate/=serviceRate;
    serviceRate=1.0;
    
}



void QueueData::calculateWardStateSpace(int wardsTotal){
    //calculates the local ward state space size
    
    Combinatorial cmb;
    vector<int> uL(wardsTotal,0);
    vector<int> lL(wardsTotal,0);
    for (int i=0; i<uL.size(); i++){
        uL[i] = capacity;
    }
    
    wardStateSpaceSize = cmb.capWithLimits(capacity,uL,lL);
    
}


void QueueData::addTimeDependency(vector<double> tDep){
    
    timeDep.resize(tDep.size());
    for (int i=0; i<timeDep.size(); i++){
        timeDep[i] = tDep[i];
    }
    
}
