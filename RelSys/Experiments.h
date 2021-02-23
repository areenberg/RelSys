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
 * File:   Experiments.h
 * Author: anders
 *
 * Created on January 26, 2021, 3:11 PM
 */

#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "WardData.h"

#include <vector>
#include <string>

class Experiments {
public:
    Experiments(int nW, WardData * wards, int sd);
    Experiments(const Experiments& orig);
    virtual ~Experiments();
    
    //METHODS
    
    //repeats the heuristic with nRep replications and
    //stores the marginal distribution of widx as a CSV-file
    void heuristicExperiment(int nRep, int widx, string label);   
    void simulationExperiment(double burnIn, double minTime, 
        int nRep, int widx, string label);   
    
    
    //VARIABLES
    vector<vector<double>> distResults; //marginal distributions (columns are replications)
    vector<long> runTimeResults; //runtime of each replication
    int memResults; //estimated maximum virtual memory consumption
    
    
    
private:

    //METHODS
    void initialize();
    void setSeed();
    double randomUniform();
    int randomInteger(int from, int to);
    vector<int> generateSeeds(int nRep);
    void saveResults(string label);
    bool numberExists(int y, vector<int> list);
    
    //VARIABLES
    int nWards, seed;
    
    mt19937 rgen; //random generator
    uniform_real_distribution<> dis;
    
    
    
    //WARD INFORMATION METHODS AND VARIABLES
    WardData * wards_pointer;
    
    //methods
    int getWardID(int ward);
    double getWardArrivalRate(int ward);
    double getWardServiceRate(int ward);
    int getWardCapacity(int ward);
    vector<double> getWardRelocationProbabilities(int ward);
    int getWardStateSpaceSize(int ward);
    void calculateWardStateSpaceSize(int ward, int numberOfWards);

    
};

#endif /* EXPERIMENTS_H */

