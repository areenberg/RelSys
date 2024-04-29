
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

