
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
#include <iomanip>
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


int main(int argc, char** argv) {
    
    
    auto start = high_resolution_clock::now(); //start time 
    
    
    //WARD SETUP
    
    //parameters
    vector<double> arrivalRates = {14.2056738,11.3679245,8.0682353,6.5424528,4.7665094,2.8983452,4.2151300,2.3915094,1.9150943,0.9433962,
                                0.9552941};
    //vector<double> arrivalRates = {11.3679245,6.5424528,4.7665094,2.8983452,2.3915094,1.9150943,0.9433962,0.9552941,26.4890391};
    vector<double> serviceRates = {0.3428191,0.2503062,0.2226145,0.7347597,0.2552630,0.2372332,0.2365980,0.1367847,0.1344521,0.1212273,0.6802304};
    //vector<double> serviceRates = {0.2503062,0.7347597,0.2552630,0.2372332,0.1367847,0.1344521,0.1212273,0.6802304,0.2773824};
    vector<int> capacity = {52,40,26,20,20,18,22,24,22,8,3};
    //vector<int> capacity = {40,20,20,18,24,22,8,3,100};
    
    //wards to include in sub-system of each test
    
    //design example
//    vector<vector<int>> incl = {{0,1,3,7,8},
//                                {0,1,8},
//                                {0,1,2,4,5,8},
//                                {0,2,3,4,5,8},
//                                {0,2,3,4,5,8},
//                                {0,2,4,5,6,8},
//                                {0,3,4,5,6,8},
//                                {0,1,2,3,7,8},
//                                {0,3,4,8}};
    
    //original setting
//    vector<vector<int>> incl = {{0,1,2,4,6},
//                                {1,2,3,5,0},
//                                {0,1,2,5,6,7},
//                                {0,1,3},
//                                {0,1,3,4,6,8},
//                                {0,1,2,5,6,7},
//                                {0,1,2,4,5,6},
//                                {0,1,2,4,6,7},
//                                {0,1,2,4,6,8},
//                                {0,1,2,5,8,9},
//                                {0,1,2,3,5,10}};
    
    //original rel. matrix
    vector<vector<double>> relProbs = {{0.0,0.08265521,0.15831102,0.008942320,0.03554432,0.14304674,0.167602077,0.18545501,0.11365490,0.022131217,0.005594279},
                                       {0.150184908,0.0,0.07687527,0.012754335,0.07673753,0.27708493,0.123951266,0.08226165,0.03398018,0.024583773,0.083098567},
                                       {0.263110997,0.05924154,0.0,0.003554453,0.01259642,0.13123587,0.203823045,0.19078884,0.09766925,0.015923550,0.003770488},
                                       {0.004814547,0.20221099,0.0,0.0,0.05339771,0.0,0.004814547,0.0,0.0,0.0,0.004814547},
                                       {0.128302181,0.07950094,0.03347470,0.0,0.0,0.04468485,0.081656734,0.07814084,0.21852669,0.002293495,0.005358570},
                                       {0.136559769,0.18739042,0.15508040,0.007864339,0.01900806,0.0,0.154181868,0.09665628,0.05443482,0.008960431,0.086519181},
                                       {0.153524196,0.06768077,0.18432780,0.010451380,0.09455126,0.16873848,0.0,0.09215592,0.11262817,0.006043909,0.005386044},
                                       {0.188539002,0.05059236,0.15901123,0.0,0.02766936,0.19394361,0.160486716,0.0,0.15485382,0.008909503,0.006066967},
                                       {0.194533623,0.09414716,0.10454087,0.0,0.06003041,0.08139556,0.106895175,0.15777673,0.0,0.171972970,0.002307912},
                                       {0.060993953,0.04898512,0.01922381,0.0,0.03844763,0.02691334,0.064079376,0.0,0.38702953,0.0,0.0},
                                       {0.012329438,0.46029902,0.03602808,0.0,0.05342756,0.02465888,0.021778166,0.0,0.0,0.0,0.0}};

    //design example
//    vector<vector<double>> relProbs = {
//        {0.00000000,0.012754335,0.07673753,0.27708493,0.08226165,0.03398018,0.024583773,0.083098567,0.351011446},
//        {0.20221099,0.000000000,0.05339771,0.00000000,0.00000000,0.00000000,0.000000000,0.004814547,0.009629095},
//        {0.07950094,0.000000000,0.00000000,0.04468485,0.07814084,0.21852669,0.002293495,0.005358570,0.243433617},
//        {0.18739042,0.007864339,0.01900806,0.00000000,0.09665628,0.05443482,0.008960431,0.086519181,0.445822041},
//        {0.05059236,0.000000000,0.02766936,0.19394361,0.00000000,0.15485382,0.008909503,0.006066967,0.508036946},
//        {0.09414716,0.000000000,0.06003041,0.08139556,0.15777673,0.00000000,0.171972970,0.002307912,0.405969663},
//        {0.04898512,0.000000000,0.03844763,0.02691334,0.00000000,0.38702953,0.000000000,0.000000000,0.144297141},
//        {0.46029902,0.000000000,0.05342756,0.02465888,0.00000000,0.00000000,0.000000000,0.000000000,0.070135682},
//        {0.07314085,0.007541375,0.03794429,0.14353755,0.17223320,0.10862249,0.017680508,0.005005639,0.000000000}};
    
    
    //for (int mWard=0; mWard<arrivalRates.size(); mWard++){
        
    //construct truncated matrix
//    vector<vector<double>> relProbs2(incl[mWard].size());
//    for (int i=0; i<incl[mWard].size(); i++){
//        relProbs2[i].resize(incl[mWard].size(),0.0);
//        for (int j=0; j<incl[mWard].size(); j++){
//            relProbs2[i][j]=relProbs[incl[mWard][i]][incl[mWard][j]];
//        }
//    }
//    
    int nWards = 11;// incl[mWard].size();
    
    WardData * wd_array = new WardData[nWards];
    
//    cout << "PARAMETERS" << endl;
//    for (int i=0; i<nWards; i++){
//        cout << "Adding ward " << incl[mWard][i] << ". Arr = " << arrivalRates[incl[mWard][i]] << ", Ser = " << 
//                serviceRates[incl[mWard][i]] << ", Cap = " << capacity[incl[mWard][i]] << endl;  
//        wd_array[i] = WardData(i,arrivalRates[incl[mWard][i]],serviceRates[incl[mWard][i]],
//                capacity[incl[mWard][i]],relProbs2[i]);
//        cout << "Relocation: ";
//        for (int j=0; j<relProbs2[i].size(); j++){
//            cout << relProbs2[i][j] << ", ";
//        }
//        cout << endl;
//    }
    for (int i=0; i<nWards; i++){
        wd_array[i] = WardData(i,arrivalRates[i],serviceRates[i],
                capacity[i],relProbs[i]);
    }
    
    //--------------------------
    //EXPERIMENTS
    //--------------------------
    
    Experiments expr(nWards,wd_array,123);
    
    //observed distribution
    vector<vector<int>> obsDists = {
                                    {0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,3,6,2,6,9,9,5,20,19,15,14,30,29,44,
                                    46,50,74,81,64,73,78,98,98,95,62,66,54,43,19},
                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,3,
                                    22,22,18,29,30,34,35,31,70,84,108,74},
                                    {0,0,0,0,0,0,0,0,0,0,0,1,2,0,2,2,1,3,4,8,17,17,43,60,74,86,90},
                                    {0,0,3,14,8,35,64,73,56,61,46,55,57,44,28,16,17,8,5,5,2},
                                    {0,0,0,0,0,0,0,2,3,9,6,21,32,26,43,43,58,38,63,54,25},
                                    {0,0,0,0,0,0,0,0,0,0,1,2,11,19,29,47,41,45,28},
                                    {0,0,0,0,0,0,0,0,0,0,0,1,0,4,4,13,18,27,41,63,48,61,51},
                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,2,17,22,35,38,31,44,19,12},
                                    {0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,4,11,17,21,35,38,36,10},
                                    {0,0,0,1,1,1,11,16,8},
                                    {27,34,17,8}
                                    };
    //expr.KSExperiment(365,1825000,50000,false,1.0,obsDists);
    int method = 0;
    expr.chiSqExperiment(365,10000,5000,false,1.0,obsDists,method);
    
    //int nRep = 1;
    
    //find widx
//    int widx = 0;
//    while (incl[mWard][widx]!=mWard){
//            widx++; 
//    }
    
    
    
    //for (int widx=0; widx<nWards; widx++){
    //    string testName = "Test" + to_string(widx+1);
    //    expr.simulationExperiment(365,1825000,nRep,widx,testName,
    //        false,1.0);
        //expr.heuristicExperiment(nRep,widx,testName);
    //}
    //}
    
    //--------------------------
    //SIMULATION
    //--------------------------
    
//    RelocSimulation * sim_pointer = new RelocSimulation[1];
//        
//    //setup model object
//    sim_pointer[0] = RelocSimulation(nWards,wd_array);
//        
//    //run calculations
//    sim_pointer[0].setSeed(1234);
//    vector<int> maxWardSamples(1,-1);
//    sim_pointer[0].simulate(365,10000,maxWardSamples,50);
//        
//    //marginal distribution
//    for (int widx=0; widx<nWards; widx++){
//        cout << "WARD " << (widx+1) << endl;
//        cout << "nWardFreDist: " << sim_pointer[0].nWardFreq[widx] << endl;
//        for (int j=0; j<sim_pointer[0].wardFreqDist[widx].size(); j++){
//            double m = sim_pointer[0].wardFreqDist[widx][j];///sim_pointer[0].nWardFreq[widx];
//            cout << m << " ";
//        }
//        cout << endl;
//    }
        
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
    
    
    return 0;
}

