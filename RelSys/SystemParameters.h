
/* 
 * File:   SystemParameters.h
 * Author: Anders Reenberg Andersen
 *
 * Created on August 11, 2021, 8:54 PM
 */

#ifndef SYSTEMPARAMETERS_H
#define SYSTEMPARAMETERS_H

#include "CustomerData.h"

#include <vector>

using namespace std;


class SystemParameters {
public:
    
    //METHODS
    double queueArrivalRate(int widx);
    double queueServiceRate(int widx);
    vector<double> queueRelProbability(int widx);
    
    //VARIABLES
    
    SystemParameters() {}; //dummy constructor
    SystemParameters(int nW, int nCusts, CustomerData * custs);
    SystemParameters(const SystemParameters& orig);
    virtual ~SystemParameters();
    
private:

    
    //METHODS
    
    
    //VARIABLES
    CustomerData * custs_pointer;
    
    //number of queues and customer types in the system
    int nQueues, nCustTypes;
    
    
};

#endif /* SYSTEMPARAMETERS_H */

