
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

class Combinatorial {
public:
    
    Combinatorial();
    Combinatorial(const Combinatorial& orig);
    virtual ~Combinatorial();
    
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






