


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

    //time-dependent arrival intensity specified as a fraction of the arrival
    //rate stated above
    vector<double> timeDep;
    
    void addTimeDependency(vector<double> tDep);
    
    //use this method to set the service rate to a value of 1
    //and scale the arrival rate keeping the load on the system
    void equalServiceRate();
    
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