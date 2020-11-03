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
 * File:   WardData.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on November 1, 2020, 7:14 PM
 */

#include "WardData.h"
#include "Combinatorial.h"

WardData::WardData(int wdnr, double arrRate, double serRate, int cap, vector<double> relProbs):
wardnr(wdnr),
arrivalRate(arrRate),
serviceRate(serRate),        
capacity(cap),
relocationProbabilities(relProbs)
{   
}

WardData::WardData(const WardData& orig) {
}

WardData::~WardData() {
}


void WardData::calculateWardStateSpace(int wardsTotal){
    //calculates the local ward state space size
    
    Combinatorial cmb;
    vector<int> uL(wardsTotal,0);
    vector<int> lL(wardsTotal,0);
    for (int i=0; i<uL.size(); i++){
        uL[i] = capacity;
    }
    
    wardStateSpaceSize = cmb.capWithLimits(capacity,uL,lL);
    
}