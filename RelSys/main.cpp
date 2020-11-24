
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

#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono; 



int main(int argc, char** argv) {
    
    auto start = high_resolution_clock::now(); //start time 
    
    int seed = 123;
    
    //--------------------------
    //WARD SETUP
    //--------------------------
    
    int nWards = 3;
    WardData * wd_array = new WardData[nWards];
    vector<vector<double>> relProbs = {{0.0,0.50,0.50},
                                       {0.50,0.0,0.50},
                                       {0.50,0.50,0.0}};
    wd_array[0] = WardData(0,1.0,0.25,5,relProbs[0]);
    wd_array[1] = WardData(1,0.75,0.40,3,relProbs[1]);
    wd_array[2] = WardData(2,1.25,0.70,3,relProbs[2]);
    
    //--------------------------
    //SIMULATION
    //--------------------------
    
    RelocSimulation * sim_pointer = new RelocSimulation[1];
    sim_pointer[0] = RelocSimulation(nWards,wd_array);

    sim_pointer[0].setSeed(seed);
    sim_pointer[0].simulate(365,365,500);
    
//    //fitting PH-parameters
//    PhaseFitter ph;
//    //parameters
//    int phases = 1;
//    int EMiterations = 100;
//    int seed = 123;
//    
//    int widx = 2; //ward index
//    ph.setInputSample(sim.blockedTimes[widx]);
//    
//    //set output PH distribution
//    ph.setHyberExponential(phases);
//    
//    //run calculations
//    ph.run(EMiterations,seed);
//    
//    //get the result
//    cout << "Distribution, pi:" << endl;
//    for (int i=0; i<phases; i++){
//        cout << ph.init_dist[i] << endl;
//    }
//    cout << "Exit-rates:" << endl;
//    for (int i=0; i<phases; i++){
//        cout << ph.exit_rate_vector[i] << endl;
//    }
    
    
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
    
    int widx, main_widx = 1;
    double arr;
    
    //create and add surrogate queues to the system
    int nhq = 2; //number of hyper queues
    HyperQueue * hq_array = new HyperQueue[nhq];
    
    int statesBlocked = 1;
    int statesOpen = 2;
    
    widx = 0;
    arr = wd_array[widx].arrivalRate*wd_array[widx].relocationProbabilities[main_widx];
    hq_array[0] = HyperQueue(widx,statesBlocked,statesOpen,
            arr,wd_array[widx].serviceRate,sim_pointer);
    hq_array[0].fitAll(seed); 
    
    widx = 2;
    arr = wd_array[widx].arrivalRate*wd_array[widx].relocationProbabilities[main_widx];
    hq_array[1] = HyperQueue(widx,statesBlocked,statesOpen,
            arr,wd_array[widx].serviceRate,sim_pointer);
    hq_array[1].fitAll(seed); 
    
    
    //specifications of main queue
    double arrivalRate = wd_array[main_widx].arrivalRate;
    double serviceRate = wd_array[main_widx].serviceRate;
    int capacity = wd_array[main_widx].capacity;
    
    //upper and lower limits in each queue. must include all queues (main and hyper queues)
    vector<int> upperLimits = {wd_array[main_widx].capacity,wd_array[main_widx].capacity,wd_array[main_widx].capacity};
    vector<int> lowerLimits = {0,0,0};
    
    HeuristicQueue q(capacity,upperLimits,lowerLimits,arrivalRate,serviceRate,nhq,hq_array);
    
    
    
    cout << "number of states = " << q.Ns << endl;
    
    //initial state distribution
    vector<double> pi(q.Ns,0); 
    srand(time(0)); double sm=0;
    for (int i=0; i<q.Ns; i++){
        pi[i] = rand()%q.Ns+1;
        sm += pi[i]; 
    }
    for (int i=0; i<q.Ns; i++){
        pi[i] /= sm;
    }
    
    //solve for steady-state distribution
    LinSolver solver;
    solver.sor(pi,q,1.0,1e-9);
    
    //print results
    cout << "MARGINAL DISTRIBUTION:" << endl;
    vector<double> x = q.marginalDist(pi);
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

