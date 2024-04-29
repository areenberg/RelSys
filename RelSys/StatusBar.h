
/* 
 * File:   StatusBar.h
 * Author: Anders Reenberg Andersen
 * 
 */

#ifndef STATUSBAR_H
#define STATUSBAR_H

#include <iostream>
#include <vector>

using namespace std;

class StatusBar {
public:
    
    StatusBar(double trgt, int sz);
    StatusBar(const StatusBar& orig);
    virtual ~StatusBar();
    
    void startBar();
    void updateBar(double val);
    void endBar();

private:

    double target, currentFrac;
    int ticksAdded;
    int ticksTotal; //length of the bar
    bool ended;

    void getCurrentFraction(double val);
};

#endif /* STATUSBAR_H */
