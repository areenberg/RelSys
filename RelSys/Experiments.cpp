/*
 * Copyright 2021 anders.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* 
 * File:   Experiments.cpp
 * Author: anders
 * 
 * Created on January 26, 2021, 3:11 PM
 */

#include "Experiments.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "WardData.h"

#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;
using namespace std::chrono;


Experiments::Experiments(int nW, WardData * wards, int sd):
wards_pointer(wards),
nWards(nW),
seed(sd)
{
    initialize();
}

Experiments::Experiments(const Experiments& orig) {
}

Experiments::~Experiments() {
}


void Experiments::initialize(){
    
    //set seed
    setSeed();
    
}



void Experiments::heuristicExperiment(int nRep, int widx, string label){
    
    cout << "Running experiment " << label << endl;
    cout << "Overall experiment seed: " << seed << endl; 
    
    //initialize vector of results
    distResults.resize(nRep);
    runTimeResults.resize(nRep,0);
    
    //seeds
    vector<int> seedList = generateSeeds(nRep);
    
    //EXPERIMENTS
    for (int rep=0; rep<seedList.size(); rep++){
        RelocEvaluation * relocEva = new RelocEvaluation[1];
        cout << "Rep " << (rep+1) << ", seed: " << seedList[rep] << endl;
        
        auto start = high_resolution_clock::now(); //start time
        
        //setup model object
        relocEva[0] = RelocEvaluation(nWards,wards_pointer);
        
        //simulate samples
        relocEva[0].runSimulation(seedList[rep],365,365,50);
    
        //evaluate ward
        relocEva[0].runHeuristic(widx);
        
        //get runtime
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        runTimeResults[rep] = duration.count();
        if (rep==0){
            memResults = relocEva[0].vmemory; //estimated maximum virtual memory usage
        }
            
        //marginal distribution
        distResults[rep].resize(relocEva[0].marginalDist.size(),0);
        for (int j=0; j<relocEva[0].marginalDist.size(); j++){
            distResults[rep][j] = relocEva[0].marginalDist[j];
        }
        delete[] relocEva; //free model from memory
        
    }
    
    //save results
    saveResults(label);
}


void Experiments::simulationExperiment(double burnIn, double minTime,
        int nRep, int widx, string label){
    
    cout << "Running experiment " << label << endl;
    cout << "Overall experiment seed: " << seed << endl; 
    
    //initialize vector of results
    distResults.resize(nRep);
    runTimeResults.resize(nRep,0);
    
    //seeds
    vector<int> seedList = generateSeeds(nRep);
    
    //EXPERIMENTS
    for (int rep=0; rep<seedList.size(); rep++){
        RelocSimulation * sim_pointer = new RelocSimulation[1];
        
        //setup model object
        sim_pointer[0] = RelocSimulation(nWards,wards_pointer);
        
        //run calculations
        sim_pointer[0].setSeed(seed);
        auto start = high_resolution_clock::now();
        
        sim_pointer[0].simulate(burnIn,minTime,50);
    
        //get runtime
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        runTimeResults[rep] = duration.count();
        if (rep==0){
            memResults = 0;
        }
        
        //marginal distribution
        distResults[rep].resize(sim_pointer[0].wardFreqDist[widx].size(),0);
        for (int j=0; j<sim_pointer[0].wardFreqDist[widx].size(); j++){
            distResults[rep][j] = sim_pointer[0].wardFreqDist[widx][j];
        }
        delete[] sim_pointer; //free model from memory
        
    }
    
    //save results
    saveResults(label);
    
}



void Experiments::setSeed(){
    mt19937 rgen(seed);
}

double Experiments::randomUniform(){
    //generate a random uniform number in the range [0,1)
    double r = dis(rgen);
    return(r);
}

int Experiments::randomInteger(int from, int to){
    int y = (int) round(randomUniform()*(to-from)+from); 
    return(y);
}

vector<int> Experiments::generateSeeds(int nRep){
    //generate nRep simulation seeds
    //for the experiments
    
    int y,from,to;
    vector<int> seedList(nRep,0);
    from = 100;
    to = 10000;
    
    
    for (int i=0; i<seedList.size(); i++){
        do {
            y = randomInteger(from,to);
        } while (numberExists(y,seedList));
        seedList[i] = y;
    }
    
    return(seedList);
}

bool Experiments::numberExists(int y, vector<int> list){
    //checks if the generated seed is already
    //on the list.
    
    bool exists = false;
    int i=0;
    while (i<list.size() && y!=list[i]){
        i++;
    }
    if (y==list[i]){
        exists = true;
    }
    
    return(exists);
}

void Experiments::saveResults(string label){
    //save experiment results in CSV-files.
    //distRes.csv contains the marginal distributions.
    //runMemRes.csv contains runtime and memory.
    
    //marginal distributions
    ofstream distResFile;
    string fName1 = "distRes_" + label + ".csv";
    distResFile.open(fName1);
    for(int i=0; i<distResults[0].size(); i++){
        for (int j=0; j<distResults.size(); j++){
            distResFile << distResults[j][i];
            if (j<(distResults.size()-1)){
                distResFile << "," << flush;
            }
        }
        distResFile << endl;
    }
    distResFile.close();
    
    //runtime and memory results
    ofstream runMemFile;
    string fName2 = "runMemRes_" + label + ".csv";
    runMemFile.open(fName2);
    for (int i=0; i<runTimeResults.size(); i++){
        runMemFile << runTimeResults[i];
        if (i<(runTimeResults.size()-1)){
            runMemFile << "," << flush;
        }
    }
    runMemFile << endl;
    runMemFile << memResults << endl;
    runMemFile.close();
    
    cout << "Results saved as " << fName1 << " and " << fName2 << endl;
    
}


