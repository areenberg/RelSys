
/* 
 * File:   Patient.cpp
 * Author: Anders Reenberg Andersen
 * 
 */

#include "Customer.h"
#include <limits>

Customer::Customer(double arr, double ser, int widx, int pidx):
arrivalClock(arr),
serviceTime(ser),
patientType(pidx),
wardTarget(widx),
active(true),
serviceClock(numeric_limits<double>::max())
{
}

Customer::Customer(const Customer& orig) {
}

Customer::~Customer() {
}
