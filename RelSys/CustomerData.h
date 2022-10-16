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
 * File:   CustomerData.h
 * Author: Anders Reenberg Andersen
 *
 * Created on August 11, 2021, 8:43 PM
 */

#ifndef CUSTOMERDATA_H
#define CUSTOMERDATA_H

#include <vector>

using namespace std;

class CustomerData {
public:
    
    //index of the queue preferred by the customer type
    int queue;
    
    //arrival rate (mean arrivals per time unit) of the customer type
    double arrivalRate;
    
    //mean length-of-stay / service time of the customer type
    double los;
    
    //vector containing the probability of relocating to an alternative queue
    //in case the preferred queue has insufficient capacity to take the customer.
    vector<double> probs;
    
    //dummy constructor (not included in cpp-file) 
    CustomerData() {};
    //parameterized constructor    
    CustomerData(int q, double arrRate, double meanLOS, vector<double> relProbs);
    CustomerData(const CustomerData& orig);
    virtual ~CustomerData();
private:

};

#endif /* CUSTOMERDATA_H */

