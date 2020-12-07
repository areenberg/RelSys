
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
#include "RelocSimulation.h"
#include "PhaseFitter.h"
#include "RelocEvaluation.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono; 



int main(int argc, char** argv) {
    
    
    auto start = high_resolution_clock::now(); //start time 
    
    //--------------------------
    //WARD SETUP
    //--------------------------
    
    //setup wards
    int nWards = 3;
    WardData * wd_array = new WardData[nWards];
    vector<vector<double>> relProbs = {{0.0,0.5,0.5},
                                       {0.5,0.0,0.5},
                                       {0.5,0.5,0.0}};
    wd_array[0] = WardData(0,1.0,0.25,3,relProbs[0]);
    wd_array[1] = WardData(1,0.75,0.40,2,relProbs[1]);
    wd_array[2] = WardData(2,1.25,0.70,3,relProbs[2]);
    
    //setup model object
    RelocEvaluation mdl(nWards,wd_array);
    
    //simulate samples
    int seed = 123;
    mdl.runSimulation(seed,365,365,50);
    
    //evaluate ward
    int widx = 2; //ward index to be evaluated
    mdl.runHeuristic(widx);
    
    cout << "MARGINAL DISTRIBUTION" << endl;
    for (int i=0; i<mdl.marginalDist.size(); i++){
        cout << mdl.marginalDist[i] << endl;
    }
    
    //--------------------------
    //EXACT SYSTEM
    //--------------------------
    
//    //<<The EntireSystem class is currently buggy!>>
//    EntireSystem sys(nWards,wd_array);
//    sys.printStateSpace();
//    
//    cout << "number of states: " << sys.nS << endl;
//    
//    //sys.selfCheck();
//    
//    vector<double> pi(sys.nS,0); 
//    srand(time(0)); double sm=0;
//    for (int i=0; i<sys.nS; i++){
//        pi[i] = rand()%sys.nS+1;
//        sm += pi[i]; 
//    }
//    for (int i=0; i<sys.nS; i++){
//        pi[i] /= sm;
//    }
//    
//    LinSolver solver;
//    solver.sorExactSystem(pi,sys,1.0,1e-9);
//    
//    
//    //PRINT RESULTS
//    cout << "MARGINAL DISTRIBUTIONS:" << endl;
//    
//    vector<vector<double>> dist(nWards);
//    for (int widx=0; widx<nWards; widx++){
//        dist[widx] = sys.marginalDist(widx,pi);
//    }
//    
//    for (int widx=0; widx<nWards; widx++){
//        for (int i=0; i<dist[widx].size(); i++){
//            cout << dist[widx][i] << endl;
//        }
//        cout << "#" << endl;
//    }
//    
//    cout << "STATS:" << endl;
//    cout << "Expected load" << endl;
//    for (int widx=0; widx<nWards; widx++){
//        cout << "Ward " << (widx+1) << " ";
//    }
//    cout << endl;
//    for (int widx=0; widx<nWards; widx++){
//        cout << sys.expectedOccupancy(widx,pi) << " ";
//    }
//    cout << endl;
//    cout << "Capacity utilization" << endl;
//    for (int widx=0; widx<nWards; widx++){
//        cout << sys.expectedOccupancyFraction(widx,pi)*100 << "%" << " ";
//    }
//    cout << endl;
    
    //--------------------------
    //HEURISTIC SYSTEM
    //--------------------------
    
//    int widx, main_widx = 1;
//    double arr;
//    
//    //create and add surrogate queues to the system
//    int nhq = 2; //number of hyper queues
//    HyperQueue * hq_array = new HyperQueue[nhq];
//    
//    int statesBlocked = 1;
//    int statesOpen = 2;
//    
//    widx = 0;
//    arr = wd_array[widx].arrivalRate*wd_array[widx].relocationProbabilities[main_widx];
//    hq_array[0] = HyperQueue(widx,statesBlocked,statesOpen,
//            arr,wd_array[widx].serviceRate,sim_pointer);
//    hq_array[0].fitAll(seed); 
//    
//    widx = 2;
//    arr = wd_array[widx].arrivalRate*wd_array[widx].relocationProbabilities[main_widx];
//    hq_array[1] = HyperQueue(widx,statesBlocked,statesOpen,
//            arr,wd_array[widx].serviceRate,sim_pointer);
//    hq_array[1].fitAll(seed); 
//    
//    
//    //specifications of main queue
//    double arrivalRate = wd_array[main_widx].arrivalRate;
//    double serviceRate = wd_array[main_widx].serviceRate;
//    int capacity = wd_array[main_widx].capacity;
//    
//    //upper and lower limits in each queue. must include all queues (main and hyper queues)
//    vector<int> upperLimits = {wd_array[main_widx].capacity,wd_array[main_widx].capacity,wd_array[main_widx].capacity};
//    vector<int> lowerLimits = {0,0,0};
//    
//    HeuristicQueue q(capacity,upperLimits,lowerLimits,arrivalRate,serviceRate,nhq,hq_array);
//    
//    
//    
//    cout << "number of states = " << q.Ns << endl;
//    
//    //initial state distribution
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
//    //solve for steady-state distribution
//    LinSolver solver;
//    solver.sor(pi,q,1.0,1e-9);
//    
//    //print results
//    cout << "MARGINAL DISTRIBUTION:" << endl;
//    vector<double> x = q.marginalDist(pi);
//    for (int i=0; i<x.size(); i++){
//        cout << x[i] << endl;
//    }
//    cout << "STATS:" << endl;
//    cout << "expected load = " << q.expectedOccupancy(pi) << endl;
//    cout << "capacity utilization = " << q.expectedOccupancyFraction(pi)*100 << "%" << endl;
//
//    auto stop = high_resolution_clock::now(); //stop time 
//    auto duration = duration_cast<milliseconds>(stop - start); 
//    cout << "Runtime: " << duration.count() << " milliseconds" << endl; 
    
    
    
    
    return 0;
}

