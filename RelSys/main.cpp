
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
    
    int nWards = 3;
    WardData * wd_array = new WardData[nWards];
    vector<vector<double>> relProbs = {{0.0,0.5,0.5},
                                       {0.5,0.0,0.5},
                                       {0.5,0.5,0.0}};
    wd_array[0] = WardData(0,1.0,0.25,5,relProbs[0]);
    wd_array[1] = WardData(1,0.75,0.40,3,relProbs[1]);
    wd_array[2] = WardData(2,1.25,0.70,3,relProbs[2]);
    
    //--------------------------
    //SIMULATION
    //--------------------------
    
    RelocSimulation sim(nWards,wd_array);
    sim.setSeed(123);
    sim.simulate(1000);
    
    //--------------------------
    //EXACT SYSTEM
    //--------------------------
    
//    EntireSystem sys(nWards,wd_array);
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
//    cout << "MARGINAL DISTRIBUTION:" << endl;
//    for (int widx=0; widx<nWards; widx++){
//        cout << "Ward " << (widx+1) << "    ";
//    }
//    cout << endl;
//    vector<vector<double>> dist(nWards);
//    for (int widx=0; widx<nWards; widx++){
//        dist[widx] = sys.marginalDist(widx,pi);
//    }
//    int capMax = 0;
//    for (int i=0; i<nWards; i++){
//        if (wd_array[i].capacity>capMax){
//            capMax = wd_array[i].capacity;
//        }
//    }
//    for (int i=0; i<=capMax; i++){
//        for (int widx=0; widx<nWards; widx++){
//            if (i<=wd_array[widx].capacity){
//                cout << dist[widx][i] << ",";
//            }else{
//                cout << " ";
//            }
//            
//        }
//        cout << endl;
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
//    //specifications of main queue
//    double arrivalRate = 1.25;
//    double serviceRate = 0.70;
//    int capacity = 3;
//    
//    //create and add surrogate queues to the system
//    int nhq = 2; //number of hyper queues
//    HyperQueue * hq_array = new HyperQueue[nhq];
//    //rates with which relocated patients arrive to the main queue during blockage
//    double arrRate[nhq] = {1.0*0.50,0.75*0.50};
//    //rates with which the relocated patients are served at the main queue
//    double serRate[nhq] = {0.25,0.40};
//    for (int i=0; i<nhq; i++){
//        int statesBlocked = 1;
//        int statesOpen = 2;
//        hq_array[i] = HyperQueue(statesBlocked,statesOpen,arrRate[i],serRate[i]);
//        //hq_array[i].fitAll(i); //fit parameters
//    }
//    hq_array[0].fitAll(0); 
//    hq_array[1].fitAll(1);
//    
//    //upper and lower limits in each queue. must include all queues (main and hyper queues)
//    vector<int> upperLimits = {3,3,3};
//    vector<int> lowerLimits = {0,0,0};
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
//    solver.sor(pi,q,1.0,1e-9);
//    
//    cout << "MARGINAL DISTRIBUTION:" << endl;
//    vector<double> x = q.marginalDist(pi);
//    for (int i=0; i<x.size(); i++){
//        cout << x[i] << endl;
//    }
//    cout << "STATS:" << endl;
//    cout << "expected load = " << q.expectedOccupancy(pi) << endl;
//    cout << "capacity utilization = " << q.expectedOccupancyFraction(pi)*100 << "%" << endl;
//
    auto stop = high_resolution_clock::now(); //stop time 
    auto duration = duration_cast<milliseconds>(stop - start); 
    cout << "Runtime: " << duration.count() << " milliseconds" << endl; 
    
    
    
    
    return 0;
}

