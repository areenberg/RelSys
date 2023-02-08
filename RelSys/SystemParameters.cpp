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
 * File:   SystemParameters.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on August 11, 2021, 8:54 PM
 */

#include "SystemParameters.h"
#include <iostream>

using namespace std;

SystemParameters::SystemParameters(int nW, int nCusts, CustomerData * custs):
custs_pointer(custs),
nCustTypes(nCusts),
nQueues(nW)
{
}

SystemParameters::SystemParameters(const SystemParameters& orig) {
}

SystemParameters::~SystemParameters() {
}


double SystemParameters::queueArrivalRate(int widx){
    //resulting arrival rate to the queue with index widx
    
    double arr=0.0;
    for (int i=0; i<nCustTypes; i++){
        if (custs_pointer[i].queue==widx){
            arr+=custs_pointer[i].arrivalRate;
        }
    }
    
    if (arr==0){
        arr = 1e-16;
    }
        
    return(arr);
}

double SystemParameters::queueServiceRate(int widx){
    //resulting service rate to the queue with index widx
    
    double totalArr = queueArrivalRate(widx);
    
    if (totalArr<=1e-16){
        return(1e16);
    }else{
        double m=0.0;
        for (int i=0; i<nCustTypes; i++){
            if (custs_pointer[i].queue==widx){
                m += custs_pointer[i].los*(custs_pointer[i].arrivalRate/totalArr);
            }
        }
        double ser = 1.0/m; //convert to rate (discharges per time unit)
        
        //cout << "Ward " << (widx+1) << " service rate: " << ser << endl;
    
        return(ser);
    }    
}

vector<double> SystemParameters::queueRelProbability(int widx){
    //resulting relocation probabilities associated with
    //the queue that has index widx
    
    double totalArr = queueArrivalRate(widx);
    vector<double> relProbs(nQueues,0);
    for (int i=0; i<nCustTypes; i++){
        if (custs_pointer[i].queue==widx){
            for (int j=0; j<relProbs.size(); j++){
                relProbs[j] += custs_pointer[i].probs[j]*(custs_pointer[i].arrivalRate/totalArr);
            }
        }
    }
    return(relProbs);
}


