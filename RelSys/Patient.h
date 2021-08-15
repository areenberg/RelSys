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


#ifndef PATIENT_H
#define PATIENT_H

using namespace std;

class Patient {
public:
    double arrivalClock, serviceTime, serviceClock;
    int patientType, wardTarget;
    //CONSTRUCTORS
    //dummy constructor (not included in cpp-file) 
    Patient() {};
    Patient(double arr, double ser, int widx, int pidx);
    Patient(const Patient& orig);
    virtual ~Patient();
    
private:
};

#endif /* PATIENT_H */