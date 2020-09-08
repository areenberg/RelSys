
/* 
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on 31. august 2020, 20.10
 */


#include "PhaseFitter.h"

#include <cstdlib>
#include <iostream>
#include<string.h>
#include <vector>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    //To do: (1) go over the remaining input density types.
    //       (2) create option to give/get input/output information as vectors.
    
    PhaseFitter ph; 
    
    
    
    //parameters
    int phases = 3;
    int EMiterations = 50;
    int seed = 1234;
    int inputDensityType = 6;
    double truncationPoint = 100;
    double dt = 1e-3; //truncationPoint/500;
    
    //input phases
    vector<vector<double>> q_in {{-2, 2, 0},
                                   {0, -2, 2},
                                   {0, 0, -2}};
    vector<double> pi_in {1, 0, 0};
    
    //set input and output types
    ph.setInputDensity(q_in,pi_in);
    ph.setSumOfExponentials(phases);
    
    ph.run(EMiterations,seed,inputDensityType,truncationPoint,dt);
    
    
    return 0;
}

