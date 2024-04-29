
/* 
 * File:   CustomerData.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on August 11, 2021, 8:43 PM
 */

#include "CustomerData.h"

CustomerData::CustomerData(int widx, double arrRate, double meanLOS, vector<double> relProbs):
queue(widx),
arrivalRate(arrRate),
los(meanLOS),
probs(relProbs)        
{
}

CustomerData::CustomerData(const CustomerData& orig) {
}

CustomerData::~CustomerData() {
}

