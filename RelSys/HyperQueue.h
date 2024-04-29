
/* 
 * File:   HyberQueue.h
 * Author: Anders Reenberg Andersen
 *
 * Created on September 15, 2020, 11:15 AM
 */

#ifndef HYBERQUEUE_H
#define HYBERQUEUE_H

#include <vector>
#include "RelocSimulation.h"


using namespace std;

class HyperQueue {
public:
    
    //variables
    vector<double> openRates; //rate at which the queue stays open
    vector<double> openDist; //probability distribution of getting one of the above rates
    vector<double> blockedRates;
    vector<double> blockedDist;
    
    //rate at which the customers from this hyper queue will be arriving to the main queue
    //when the hyper queue is blocked, i.e. arrivalRate = lambda_hq * p_{hyperqueue, mainqueue}
    double arrivalRate;
    //rate at which the customers are discharged in the main queue.
    double serviceRate;
    
    int numberOfStates;
    
    //methods
    void fitOpenPH(int seed);
    void fitBlockedPH(int seed);
    void fitAll(int seed);
  
    //dummy constructor (not included in cpp-file) 
    HyperQueue() {};
    //parameterized constructor
    HyperQueue(int widx, int statesBlocked, int statesOpen, double aRate, double sRate, RelocSimulation *sm);
    
    HyperQueue(const HyperQueue& orig);
    virtual ~HyperQueue();
    
    
private:

    RelocSimulation *sim_pointer;
    int wardIndex;
  
    void lowLoadOpenPH();
    void lowLoadBlockedPH();
    
};

#endif /* HYBERQUEUE_H */

