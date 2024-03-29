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
