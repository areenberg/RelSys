

#ifndef CUSTOMER_H
#define CUSTOMER_H

using namespace std;

class Customer {
public:
    double arrivalClock, serviceTime, serviceClock;
    int patientType, wardTarget;
    bool active;
    //CONSTRUCTORS
    //dummy constructor (not included in cpp-file) 
    Customer() {};
    Customer(double arr, double ser, int widx, int pidx);
    Customer(const Customer& orig);
    virtual ~Customer();
    
private:
};

#endif /* CUSTOMER_H */