
/* 
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on September 14, 2020, 7:38 PM
 */

#include "HeuristicQueue.h"
#include "HyperQueue.h"
#include "LinSolver.h"
#include "WardData.h"
#include "EntireSystem.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono; 



int main(int argc, char** argv) {
    
    auto start = high_resolution_clock::now(); //start time 
    
    //EXACT SYSTEM
    int nWards = 2;
    WardData * wd_array = new WardData[nWards];
    vector<vector<double>> relProbs = {{0.0,1.0},
                                       {1.0,0.0}};
    wd_array[0] = WardData(0,2.0,0.5,5,relProbs[0]);
    wd_array[1] = WardData(1,1.5,0.2,10,relProbs[1]);
    
    EntireSystem sys(nWards,wd_array);
    
    
      //HEURISTIC
//    //specifications of main queue
//    double arrivalRate = 2.890697;
//    double serviceRate = 0.4447632;
//    int capacity = 15;
//    
//    //create and add surrogate queues to the system
//    int nhq = 4; //number of hyper queues
//    HyperQueue * hq_array = new HyperQueue[nhq];
//    //rates with which relocated patients arrive to the main queue during blockage
//    double arrRate[nhq] = {3.585529*0.25,3.709466*0.25,2.546503*0.25,3.605887*0.25};
//    //rates with which the relocated patients are served at the main queue
//    double serRate[nhq] = {0.1364088,0.6260197,0.4431681,0.3161356};
//    for (int i=0; i<nhq; i++){
//        int statesBlocked = 1;
//        int statesOpen = 2;
//        hq_array[i] = HyperQueue(statesBlocked,statesOpen,arrRate[i],serRate[i]);
//        hq_array[i].fitAll(i); //fit parameters
//    }
//    
//    //upper and lower limits in each queue. must include all queues (main and hyper queues)
//    vector<int> upperLimits = {15,10,10,10,10};
//    vector<int> lowerLimits = {0,0,0,0,0};
//    
//    HeuristicQueue q(capacity,upperLimits,lowerLimits,arrivalRate,serviceRate,nhq,hq_array);
//    
//    
//    
//    cout << "number of states = " << q.Ns << endl;
//    
//    vector<double> pi(q.Ns,0); 
//    srand(time(0)); double sm=0;
//    for (int i=0; i<q.Ns; i++){
//        pi[i] = rand()%q.Ns+1;
//        sm += pi[i]; 
//    }
//    for (int i=0; i<q.Ns; i++){
//        pi[i] /= sm;
//    }
//    
//    LinSolver solver;
//    solver.sor(pi,q,1.0,1e-6);
//    //solver.sorOnDemand(pi,q,1.0,1e-6);
//    //solver.powerMethod(pi,q,1e-6);
//    
//    cout << "MARGINAL DISTRIBUTION:" << endl;
//    vector<double> x = q.marginalDist(pi);
//    double occ = 0;
//    for (int i=0; i<x.size(); i++){
//        cout << x[i] << endl;
//    }
//    cout << "STATS:" << endl;
//    cout << "expected load = " << q.expectedOccupancy(pi) << endl;
//    cout << "capacity utilization = " << q.expectedOccupancyFraction(pi)*100 << "%" << endl;
//    cout << "rejection probability = " << q.rejectionProbability(pi) << endl; 
//    
    auto stop = high_resolution_clock::now(); //stop time 
    auto duration = duration_cast<milliseconds>(stop - start); 
    cout << "Runtime: " << duration.count() << " milliseconds" << endl; 
    
    
    
    
    return 0;
}

