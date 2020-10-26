
/* 
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on September 14, 2020, 7:38 PM
 */

#include "Queue.h"
#include "HyperQueue.h"
#include "LinSolver.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono; 


/*
 * 
 */


void normalize(vector<double> &pi){
    
    double sm = 0;
    for (int i=0; i<pi.size(); i++){
        sm += pi[i];
    }
    for (int i=0; i<pi.size(); i++){
        pi[i] /= sm;
    }
    
}



int main(int argc, char** argv) {
    
    //specifications of main queue
    double arrivalRate = 2.0;
    double serviceRate = 0.5;
    int capacity = 20;
    
    //create and add surrogate queues to the system
    int nhq = 2; //number of hyper queues
    HyperQueue * hq_array = new HyperQueue[nhq];
    //rates with which relocated patients arrive to the main queue during blockage
    double arrRate[nhq] = {4.5*0.5,1.5*0.5};
    //rates with which the relocated patients are served at the main queue
    double serRate[nhq] = {0.50,0.20};
    for (int i=0; i<nhq; i++){
        int statesBlocked = 1;
        int statesOpen = 2;
        hq_array[i] = HyperQueue(statesBlocked,statesOpen,arrRate[i],serRate[i]);
        hq_array[i].fitAll(i); //fit parameters
    }
    
    vector<int> upperLimits = {20,20,20};
    vector<int> lowerLimits = {0,0,0};
    
    Queue q(capacity,upperLimits,lowerLimits,arrivalRate,serviceRate,nhq,hq_array);
    
    
    auto start = high_resolution_clock::now(); //start time 
    
    cout << "number of states = " << q.Ns << endl;
    
    
//    for (int i=0; i<803950; i++){
//        //q.allIngoing();
//        q.nextCurrentState();
//    }
//    
//    q.allIngoing();
//    for (int j=0; j<q.fromIdxSize; j++){
//        cout << q.jumpFromIdx[j] << " "; 
//    }
    //notes: in sidx=803950, jump from jidx=804922 causes an error due to out of bound.
    //find out how this jump was derived.
    
    vector<double> pi(q.Ns,0); 
    for (int i=0; i<q.Ns; i++){
        pi[i] = 1.0/q.Ns;
    }
    
    LinSolver solver;
    solver.sor(pi,q,1.0,1e-9);
//    solver.sorOnDemand(pi,q,1.0,1e-6);
//    solver.powerMethod(pi,q,1e-6);
    
    cout << "MARGINAL DISTRIBUTION:" << endl;
    vector<double> x = q.marginalDist(pi);
    double occ = 0;
    for (int i=0; i<x.size(); i++){
        cout << x[i] << endl;
    }
    cout << "STATS:" << endl;
    cout << "expected load = " << q.expectedOccupancy(pi) << endl;
    cout << "capacity utilization = " << q.expectedOccupancyFraction(pi)*100 << "%" << endl;
    
    auto stop = high_resolution_clock::now(); //stop time 
    auto duration = duration_cast<milliseconds>(stop - start); 
    cout << "Runtime: " << duration.count() << " milliseconds" << endl; 
    
    
    
    
    return 0;
}

