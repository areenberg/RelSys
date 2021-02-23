
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
#include "Experiments.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono; 


vector<vector<double>> uniformRelocMatrix(int nWards){
    //generates a relocation matrix with uniform
    //patient relocation
    
    cout << "Relocation matrix:" << endl;
    double relprob = 1.0/((double)nWards-1.0);
    vector<vector<double>> mat(nWards);
    for (int i=0; i<nWards; i++){
        mat[i].resize(nWards,0.0);
        for (int j=0; j<mat[i].size(); j++){
            if (i!=j){
                mat[i][j] = relprob;
            }
            cout << mat[i][j] << " " << flush;
        }
        cout << endl;
    }
    
    return(mat);
}


vector<vector<double>> neighborRelocMatrix(int nWards){
    //generates a relocation matrix where patients
    //are uniformly distributed between the two
    //neighboring wards.
    //assumes nWards>2.
    
    cout << "Relocation matrix:" << endl;
    double relprob = 0.5;
    int idx;
    vector<vector<double>> mat(nWards);
    for (int i=0; i<nWards; i++){
        mat[i].resize(nWards,0.0);
        for (int j=1; j<=2; j++){
            idx = i+j;
            if (idx>=nWards){
                idx -= nWards;    
            }
            mat[i][idx] = relprob;
        }
    }
    
    
    for (int i=0; i<nWards; i++){
        for (int j=0; j<nWards; j++){
            cout << mat[i][j] << " " << flush;
        }
        cout << endl;
    }
    
    return(mat);
}

int main(int argc, char** argv) {
    
    
    auto start = high_resolution_clock::now(); //start time 
    
    int beds = 15;
    double arr = 0.80 * (0.1*beds);
    
    //WARD SETUP
    int nWards = 6;
    WardData * wd_array = new WardData[nWards];
    vector<vector<double>> relProbs = uniformRelocMatrix(nWards);
    
    wd_array[0] = WardData(0,arr,0.1,beds,relProbs[0]);
    wd_array[1] = WardData(1,arr,0.1,beds,relProbs[1]);
    wd_array[2] = WardData(2,arr,0.1,beds,relProbs[2]);
    wd_array[3] = WardData(3,arr,0.1,beds,relProbs[3]);
    wd_array[4] = WardData(4,arr,0.1,beds,relProbs[4]);
    wd_array[5] = WardData(5,arr,0.1,beds,relProbs[5]);
    
    //--------------------------
    //EXPERIMENTS
    //--------------------------
    
    Experiments expr(nWards,wd_array,123);
    
    int nRep = 1;
    int widx = 0;
    double burnIn = 365;
    double minTime = 365;
    expr.simulationExperiment(burnIn,minTime,nRep,widx,"test");
    //expr.heuristicExperiment(nRep,widx,"Test15");
    
    //--------------------------
    //HEURISTIC
    //--------------------------
    
    //setup model object
//    RelocEvaluation mdl(nWards,wd_array);
//    
//    //simulate samples
//    int seed = 123;
//    mdl.runSimulation(seed,365,365,50);
//    
//    //evaluate ward
//    int widx = 0; //ward index to be evaluated
//    mdl.runHeuristic(widx);
//    
//    cout << "MARGINAL DISTRIBUTION" << endl;
//    for (int i=0; i<mdl.marginalDist.size(); i++){
//        cout << mdl.marginalDist[i] << endl;
//    }
//    
//    cout << "Memory consumption: " << mdl.vmemory << "KB" << endl;
    
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
    
    
    
    return 0;
}

