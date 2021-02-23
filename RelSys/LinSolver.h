/*
 * Copyright 2020 Anders Reenberg Andersen.
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
 * File:   LinSolver.h
 * Author: Anders Reenberg Andersen
 * 
 * Created on September 24, 2020, 1:05 PM
 */


#ifndef LINSOLVER_H
#define LINSOLVER_H

#include "HeuristicQueue.h"
#include "EntireSystem.h"

#include <vector>
#include <iostream>


using namespace std;

class LinSolver {
public:
    
    void sor(vector<double> &pi, HeuristicQueue &q, double relaxation, double eps); //derive the state distribution using SOR. requires the transition matrix is stored.
    void sorOnDemand(vector<double> &pi, HeuristicQueue &q, double relaxation, double eps); //employ SOR using on-demand calculations of the state transitions
    void powerMethod(vector<double> &pi, HeuristicQueue &q, double eps); //employs the power method. requires the transition matrix is stored. 
    
    //SOR algorithm for the exact system
    void sorExactSystem(vector<double> &pi, EntireSystem &q, double relaxation, double eps);
    
    double vmemory;
    
    LinSolver();
    LinSolver(const LinSolver& orig);
    virtual ~LinSolver();
    
private:

    void scale(vector<vector<double>> &values);
    void normalize(vector<double> &pi);
    void embeddedChain(vector<vector<double>> &values);
    double memUsage();
    
};

#endif /* LINSOLVER_H */
