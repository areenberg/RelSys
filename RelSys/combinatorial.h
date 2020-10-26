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
 * File:   combinatorial.h
 * Author: Anders Reenberg Andersen
 *
 * Created on 7. marts 2020, 23:41
 */

#ifndef COMBINATORIAL_H
#define COMBINATORIAL_H

#include <vector>

using namespace std;

class combinatorial {
public:
    
    combinatorial();
    combinatorial(const combinatorial& orig);
    virtual ~combinatorial();
    
    //permutation methods
    int capWithLimits(int &K, vector<int> &upperLim, vector<int> &lowerLim);
    //int capWeights(int capacity, int nodes);
    //int capWithUpperBound(vector<int> &limits);
    //int capWithoutWeights_1(int capacity, int nodes);
    //int capWithoutWeights_2(int capacity, int nodes); //faster but only stable for few nodes
    //long double numberAmongOwners(int amount, int owners); //only stable for few owners
   
private:

    //basic functions
    //long double factorial(int x);
    //void printJ(int lng, int j[]); //used in tests
    
    //parameters
    int m, K_use;
    

};

#endif /* COMBINATORIAL_H */






