/*
 * Copyright 2021 Anders Reenberg Andersen.
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

