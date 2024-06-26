
/* 
 * File:   LinSolver.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on September 24, 2020, 1:05 PM
 */

#include "HeuristicQueue.h"
#include "LinSolver.h"
#include "StatusBar.h"

#include <vector>
//#include <cstdlib>
#include <cmath>
//#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <random>
#include <limits.h>

LinSolver::LinSolver() {
    
    
    
    //set maximum generated number
    if (RAND_MAX==numeric_limits<int>::max()){
        mxRnd = RAND_MAX;
    }else{
        mxRnd = RAND_MAX+1;
    }

}

LinSolver::LinSolver(const LinSolver& orig) {
}

LinSolver::~LinSolver() {
}


void LinSolver::sor(vector<double> &pi, HeuristicQueue * q, double relaxation, double eps){
    
    cout << "Engaging SOR." << endl;
    
    double sm, x_new, diff, tol;
    int maxIter = 1e5, iter = 0;
    
    cout << "Building transition matrix..." << endl;
    cout << "Step 1:" << endl;
    q->buildTransposedChain(); //generate and store the entire transposed transition matrix
    
    cout << "Step 2:" << endl;
    scale(q->qValues); //scale the transposed transition matrix
    
    //get memory consumption at runtime (for Linux operating systems)
    vmemory = memUsage();
    
    cout << "Solving state distribution..." << endl;
    StatusBar sbar(1e-4-eps,30);
    do {
        tol = 0.0;
        for (int i=0; i<q->Ns; i++){
            sm = 0;
            for (int j=0; j<q->qColumnIndices[i].size(); j++){
                sm += q->qValues[i][j] * pi[q->qColumnIndices[i][j]];
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
            //cout << tol << endl;
            double bval = 1e-4-tol;
            sbar.updateBar(bval);
        }
        
    }while (tol>eps && iter<maxIter);
    sbar.endBar();
    cout << "SOR stats: " << iter << " iterations, tolerance " << tol << endl;
    
}

//void LinSolver::sorOnDemand(vector<double> &pi, HeuristicQueue &q, double relaxation, double eps){
//    
//    cout << "Engaging SOR (on demand)." << endl;
//    
//    double sm, x_new, diff, tol;
//    int maxIter = 1e4, iter = 0;
//    
//    do {
//        tol = 0.0;
//        q.initializeState();
//        for (int i=0; i<q.Ns; i++){
//            //cout << "state " << i << endl;
//            q.allIngoing();
//            //scale
//            double scaler = 1.0/q.jumpFromRate[q.fromIdxSize-1];
//            for (int j=0; j<q.fromIdxSize; j++){
//                q.jumpFromRate[j] *= scaler; 
//            }
//            
//            sm = 0;
//            for (int j=0; j<q.fromIdxSize; j++){
//                sm += q.jumpFromRate[j] * pi[q.jumpFromIdx[j]];
//            }
//            
//            x_new = pi[i] - relaxation*sm;
//            
//            diff = abs(x_new-pi[i])/pi[i];
//                
//            if (diff>tol){
//                tol = diff;
//            }
//            pi[i] = x_new;
//            q.nextCurrentState();
//        }
//        
//        normalize(pi);
//        
//        iter++;
//        cout << tol << endl;
//    }while (tol>eps && iter<maxIter);
//    cout << "SOR stats: " << iter << " iterations, tolerance " << tol << endl;
//    
//}


void LinSolver::scale(vector<vector<double>> &values){
    //scales the transposed and stored transition matrix so all diagonal elements equal 1
    
    double scaler;
    int scalerIdx;
    
    StatusBar sbar(values.size(),30);
    for (int i=0; i<values.size(); i++){
        scalerIdx=0;
        while(values[i][scalerIdx]>=0){
            scalerIdx++;
        }
        scaler = 1.0/values[i][scalerIdx];
        for (int j=0; j<values[i].size(); j++){
            values[i][j] *= scaler;
        }
        if (i%100==0){
            double bval = i;
            sbar.updateBar(bval);
        }
    }
    sbar.endBar();
    
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

void LinSolver::powerMethod(vector<double> &pi, HeuristicQueue * q, double eps){
    
    cout << "Engaging power method." << endl;
    
    double x_new, diff, tol;
    int maxIter = 1e4, iter = 0;
    
    q->buildTransposedChain(); //generate and store the entire transposed transition matrix
    
    embeddedChain(q->qValues); //translate into the embedded chain
    
    do {
        tol = 0.0;
        for (int i=0; i<q->Ns; i++){
            
            x_new = 0;
            for (int j=0; j<q->qColumnIndices[i].size(); j++){
                x_new += q->qValues[i][j] * pi[q->qColumnIndices[i][j]];
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
        if (iter%10==0){
            cout << tol << endl;
        }
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

double LinSolver::memUsage(){
    //memory usage in kilobytes at runtime for Linux operating systems
    //from tutorialspoint.com, How to get memory usage at runtime using C++? 
    
   double vmem = 0.0;
   //resident_set = 0.0;
   ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc directory
   //create some variables to get info
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
   stat_stream.close();
   //long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured
   //to use 2MB pages
   vmem = vsize / 1024.0;
   //resident_set = rss * page_size_kb;
    
   return(vmem);
}

double LinSolver::randomUniform(){
    //generate a random uniform number in the range [0,1)
    return((double)rand()/((double)mxRnd));
}

double LinSolver::randomExponential(double rate){
    //generate a random exponential double
    return(log(1-randomUniform())/(-rate));
}
