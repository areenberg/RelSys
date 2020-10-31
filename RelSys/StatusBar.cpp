/*
 * Copyright 2020 anders.
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
 * File:   StatusBar.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on October 28, 2020, 12:59 PM
 */

#include "StatusBar.h"

#include <iostream>
#include <cmath>

using namespace std;


StatusBar::StatusBar(double trgt, int sz):
    target(trgt),
    ticksTotal(sz),
    ticksAdded(0),
    currentFrac(0.0),
    ended(false)
{
    startBar();
}

StatusBar::StatusBar(const StatusBar& orig) {
}

StatusBar::~StatusBar() {
}


void StatusBar::startBar(){
    
    cout << "||" << flush;
    
}

void StatusBar::updateBar(double val){
    
    getCurrentFraction(val);
    
    int ticksRequired = floor(currentFrac*ticksTotal);
    if (ticksRequired>ticksAdded){
        int diff = ticksRequired-ticksAdded;
        for (int i=0; i<diff; i++){
            cout << "#" << flush;
            
        }
        ticksAdded = ticksRequired;
    }
    
}

void StatusBar::endBar(){
    
    if (ended==false){
        if (currentFrac<1.0 && ticksAdded<ticksTotal){
            currentFrac = 1.0;
            int diff = ticksTotal-ticksAdded;
            for (int i=0; i<diff; i++){
                cout << "#";
            }
        }
        cout << "||" << endl;
        ended = true;
    }
    
}

void StatusBar::getCurrentFraction(double val){
    
    currentFrac = val/target;
    
}



