
/* 
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on 31. august 2020, 20.10
 */

#include <cstdlib>
#include <iostream>
#include<string.h>

#include "PhaseFitter.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    //To do: (1) go over the remaining input density types.
    //       (2) create option to give/get input/output information as vectors.
    
    PhaseFitter ph; 
    
    //run fitting program
    int phases = 3;
    int EMiterations = 50;
    int seed = 1234;
    int sampleOrDensity = 2;
    int inputDensityType = 6;
    double truncationPoint = 100;
    double dt = 1e-3; //truncationPoint/500;
    
    ph.setSumOfExponentials(phases);
    
    ph.run(EMiterations,seed,sampleOrDensity,inputDensityType,truncationPoint,dt);
    
    
    return 0;
}

