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



#ifndef QUEUEDATA_H
#define QUEUEDATA_H

#include <iostream>
#include <vector>

using namespace std;

class QueueData {
public:
    
    //nr and capacity of the ward
    //nr corresponds to the index in the ward array
    int wardnr, capacity;

    int wardStateSpaceSize; //state space size for this single ward
    
    //characteristics of patients associated with the ward
    double arrivalRate, serviceRate;

    //fraction of items/customers/patients that are relocated to other wards during shortage
    //this vector corresponds to an entire row in the ward-patients array, for instance:
    //for three wards with indices 0,1,2 and uniform relocation of patients, ward 1
    //would have relocationProbabilities = {0.5,0.0,0.5};
    vector<double> relocationProbabilities; 

    //calculate local state space size
    void calculateWardStateSpace(int wardsTotal);
    
    //dummy constructor (not included in cpp-file) 
    QueueData() {};
    //parameterized constructor
    QueueData(int wdnr, double arrRate, double serRate, int cap, vector<double> relProbs);
    QueueData(const QueueData& orig);
    virtual ~QueueData();
    
private:

};

#endif /* WARDDATA_H */