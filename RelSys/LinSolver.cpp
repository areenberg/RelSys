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
 * File:   LinSolver.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on September 24, 2020, 1:05 PM
 */

#include "Queue.h"
#include "LinSolver.h"

#include <iostream>
#include <vector>
#include <cstdlib>

LinSolver::LinSolver() {
    
    
    
}

LinSolver::LinSolver(const LinSolver& orig) {
}

LinSolver::~LinSolver() {
}


void LinSolver::sor(vector<double> &pi, Queue &q, double relaxation, double eps){
    
    double sm, x_new, diff, tol;
    int maxIter = 1e4, iter = 0;
    
    q.buildTransposedChain(); //generate and store the entire transposed transition matrix
    
    scale(q.qValues); //scale the transposed transition matrix
    
    do {
        tol = 0.0;
        for (int i=0; i<q.Ns; i++){
            sm = 0;
            for (int j=0; j<q.qColumnIndices[i].size(); j++){
                sm += q.qValues[i][j] * pi[q.qColumnIndices[i][j]];
            }
            
            x_new = pi[i] - relaxation*sm;
            diff = abs(x_new-pi[i])/pi[i];
            if (diff>tol){
                tol = diff;
            }
            pi[i] = x_new;
        }
        
        normalize(pi);
        
        iter++;
        if (iter%10==0){
            cout << tol << endl;
        }
        
    }while (tol>eps && iter<maxIter);
    cout << "SOR stats: " << iter << " iterations, tolerance " << tol << endl;
    
}

void LinSolver::sorOnDemand(vector<double> &pi, Queue &q, double relaxation, double eps){
    
    double sm, x_new, diff, tol;
    int maxIter = 1e4, iter = 0;
    
    do {
        tol = 0.0;
        q.initializeState();
        for (int i=0; i<q.Ns; i++){
            //cout << "state " << i << endl;
            q.allIngoing();
            //scale
            double scaler = 1.0/q.jumpFromRate[q.fromIdxSize-1];
            for (int j=0; j<q.fromIdxSize; j++){
                q.jumpFromRate[j] *= scaler; 
            }
            
            sm = 0;
            for (int j=0; j<q.fromIdxSize; j++){
                sm += q.jumpFromRate[j] * pi[q.jumpFromIdx[j]];
            }
            
            x_new = pi[i] - relaxation*sm;
            
            diff = abs(x_new-pi[i])/pi[i];
                
            if (diff>tol){
                tol = diff;
            }
            pi[i] = x_new;
            q.nextCurrentState();
        }
        
        normalize(pi);
        
        iter++;
        cout << tol << endl;
    }while (tol>eps && iter<maxIter);
    cout << "SOR stats: " << iter << " iterations, tolerance " << tol << endl;
    
}


void LinSolver::scale(vector<vector<double>> &values){
    //scales the transposed and stored transition matrix so all diagonal elements equal 1
    
    double scaler;
    
    for (int i=0; i<values.size(); i++){
        scaler = 1.0/values[i][values[i].size()-1];
        for (int j=0; j<values[i].size(); j++){
            values[i][j] *= scaler;
        }
    }
    
}


void LinSolver::normalize(vector<double> &pi){
    
    double sm = 0;
    for (int i=0; i<pi.size(); i++){
        sm += pi[i];
    }
    for (int i=0; i<pi.size(); i++){
        pi[i] /= sm;
    }
    
}

void LinSolver::powerMethod(vector<double> &pi, Queue &q, double eps){
    
    double x_new, diff, tol;
    int maxIter = 1e4, iter = 0;
    
    q.buildTransposedChain(); //generate and store the entire transposed transition matrix
    
    embeddedChain(q.qValues); //translate into the embedded chain
    
    do {
        tol = 0.0;
        for (int i=0; i<q.Ns; i++){
            
            x_new = 0;
            for (int j=0; j<q.qColumnIndices[i].size(); j++){
                x_new += q.qValues[i][j] * pi[q.qColumnIndices[i][j]];
            }
            
            //for convergence validation
            diff = abs(x_new-pi[i])/pi[i];
            if (diff>tol){
                tol = diff;
            }
            pi[i] = x_new;
        }
        
        normalize(pi);
        
        iter++;
    }while (tol>eps && iter<maxIter);
    cout << "PM stats: " << iter << " iterations, tolerance " << tol << endl;
    
}

void LinSolver::embeddedChain(vector<vector<double>> &values){
    //translates the stored transition matrix into the embedded Markov chain
    
    double delta_t = 0;
    for (int i=0; i<values.size(); i++){ //derive the largest absolute diagonal element
        if (values[i][values[i].size()-1]<delta_t){
            delta_t = values[i][values[i].size()-1];
        }
    }
    delta_t = (1 / abs(delta_t))-1e-3;
    for (int i=0; i<values.size(); i++){ //scale all transitions
        for (int j=0; j<values[i].size(); j++){
            values[i][j] *= delta_t;
        }
    }
    for (int i=0; i<values.size(); i++){ //add 1 to diagonal
        values[i][values[i].size()-1] += 1.0;
    }
    
}