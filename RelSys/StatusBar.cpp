
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



