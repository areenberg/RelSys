# RelSys (**Rel**ocation **Sys**tem)
RelSys is a tool for evaluating a system of queues where arriving customers can be relocated to an alternative queue if none of the servers in the preferred queue are idle.
The source code is written in C++, but we have developed a module for the users preferring Python (Linux).

# Table of contents

1. Description of RelSys
   * Input parameters
   * Output types
2. How to use
   * Python (Linux)
   * C++
3. How to cite
4. Licence

# Description of RelSys

Consider a number of parallel queues where the capacity of the queue equals the number of servers. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue (also denoted the Erlang loss or Erlang-B system), a customer is rejected and lost from the system if all the servers are occupied upon arrival. In RelSys, we allow customers to be transferred with a probability to one of the other queues. If there is an idle server in the alternative queue, the customer is accepted and served with an exponentially distributed time with the same rate-parameter as the customer would have in the preferred queue. The figure below depicts an example featuring two queues where customers are relocated (i.e. transferred) with a probability to the other queue whenever the preferred queue is full. 

<img src="https://github.com/areenberg/RelSys/blob/development/images/example_system.jpeg" width="430" height="500">

## Input parameters

RelSys has five types of input parameters for the model:

* An arrival rate vector. Each element corresponds to a customer type.
* A service time vector. Each element corresponds to a customer type.
* A relocation probability matrix. Rows correspond to customer types and columns to queues.  
* A capacity vector. Each element corresponds to a queue.
* A preferrence vector. Each element indicates the preferred queue of each customer type.

## Output types

RelSys has six output types:

* Occupancy *probability* distributions.
* Occupancy *frequency* distributions (only available if the system was evaluated using simulation).
* Shortage probabilities.
* Availability probabilities.
* Expected occupancy.
* Expected fraction of capacity occupied.

# How to use

RelSys is available for Python users on Linux, and for C++ users on all operating systems. 

## Python (Linux)

We have created a RelSys module for Python with `pybind11`. Head to the directory `Python/Linux/`, or run `wget https://github.com/areenberg/RelSys/blob/development/Python/Linux/relsys.cpython-310-x86_64-linux-gnu.so` to download the SO-file for the module.

Start by importing the module,

```python
import relsys
```

Now, specify the input parameters for the model. In this example, we consider a system containing 4 customer types and 4 queues.

```python
arrivalRates = [0.8,2.5,0.6,2.8]

serviceTimes = [10,5,10,8]

capacity = [15,20,10,30]

relocationProbabilities = [[0.0,0.4,0.1,0.5],
                           [0.3,0.0,0.5,0.0],
                           [0.0,0.5,0.0,0.5],
                           [0.2,0.3,0.5,0.0]]

preferredQueue = [0,1,2,3]
```

The input parameters are imported using the `input` function,

```python
relsys.input(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)
```

The model can now be evaluated with `run`,

```python
relsys.run()
```

Return the resulting occupancy distributions with `getDensity` and shortage probabilities with `getShortageProb`,

```python
for queueIdx in range(4):
    print(relsys.getDensity(queueIdx))

for queueIdx in range(4):
    print(relsys.getShortageProb(queueIdx))
```

### The complete example

```python
#import the module
import relsys 

#arrival rates of each customer type
arrivalRates = [0.8,2.5,0.6,2.8]

#mean service time of each customer type
serviceTimes = [10,5,10,8]

#capacity of each queue
capacity = [15,20,10,30]

#fraction of rejected customers that are moved to an alternative queue node
#this is a number of customers x number of queues matrix
relocationProbabilities = [[0.0,0.4,0.1,0.5],
                           [0.3,0.0,0.5,0.0],
                           [0.0,0.5,0.0,0.5],
                           [0.2,0.3,0.5,0.0]]

#queue indices preferred by each customer type
preferredQueue = [0,1,2,3]

#import the parameters
relsys.input(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)

#run the model
relsys.run()

#check the resulting occupancy distribution of each queue 
for queueIdx in range(4):
    print(relsys.getDensity(queueIdx))

#check the resulting shortage probabilities of each queue 
for queueIdx in range(4):
    print(relsys.getShortageProb(queueIdx))
```

## List of functions

### Import data
* `input`. Import the input data to the model.

### Model settings

* `setType`. Set the method to use in the evaluation of the model ("simulation" (default), "approximation", "auto"). 
* `queuesEval`. Set the indices of queues to evaluate.
* `equalize`. Specify if service times should be equalized and loads correspondingly adjusted (True: On, False: Off (default)).
* `setVerbose`. Control verbose (True: On, False: Off (default)).
* `setSeed`. Set the seed.
* `setBurnIn`. Set the burn-in time of the simulation.
* `setSimTime`. Set the simulation time.
* `setSamples`. Set the minimum number of open/shortage samples.
* `setHyperPhases`. Set the number of phases in the hyper-exponential distributions accounting for the open/shortage time.

### Run calculations

* `run`. Evaluate the model using the input parameters.

### Get results

* `getDensity`. Return the density distribution of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getFreq`. Return the frequency distribution of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getShortageProb`. Return the shortage probability of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getAvailProb`. Return the probability that at least one server is available. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getExpOccupany`. Return the expected number of occupied servers.
* `getExpOccFraction`. Return the expected fraction of occupied servers.

### Return imported variables

* `getArrivalRates`. Return the imported arrival rates.
* `getServiceTimes`. Return the imported service times.
* `getCapacity`. Return the imported capacities.
* `getReloc`. Return the imported relocation probabilities.
* `getPreferredQueue`. Return the imported preferred queues.

## C++

The approach in C++ is very similar to that of Python. The directory `RelSys/` contains the complete source code for RelSys. Start by heading here. Create a `main.cpp` file (or modify the example `main.cpp` that is already in the directory).

Write the following into the `main.cpp` file,

```c++
#include "Model.h"

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

    //model code goes here

    return 0;
}
```

Now, put in the parameters of the model. Once again, we consider a system containing 4 customer types and 4 queues. Note that in C++, we also have to specify the indices of the queues we want to evaluate,

```c++
    //arrival rates for each customer type
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};
    
    //mean service time for each customer type
    vector<double> serviceTimes = {10,5,10,8};

    //capacity of each queue
    vector<int> capacity = {15,20,10,30};
    
    //fraction of rejected customers that are moved to an alternative queue node
    //this is an nCustomerTypes x nQueues matrix
    vector<vector<double>> relocationProbabilities = {{0.0,0.4,0.1,0.5},
                                                      {0.3,0.0,0.5,0.0},
                                                      {0.0,0.5,0.0,0.5},
                                                      {0.2,0.3,0.5,0.0}};
    
    //queue indices preferred by each customer type
    vector<int> preferredQueue = {0,1,2,3};

    //indices of the queues to be evaluated    
    vector<int> evaluatedQueue = {0,1,2,3};

```

Now we can create and evaluate the model,

```c++
   //create the model object
    Model mdl(arrivalRates,serviceTimes,capacity,
            relocationProbabilities,preferredQueue,
            evaluatedQueue);
    
    //now evaluate the model
    mdl.runModel();
```

The following returns the resulting occupancy distributions and shortage probabilities,

```c++
    cout << "Occupancy distributions:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << "----" << "Queue " << evaluatedQueue[i] << "----" << endl;
        for (int j=0; j<mdl.queueDenDist[evaluatedQueue[i]].size(); j++){
            cout << mdl.queueDenDist[evaluatedQueue[i]][j] << endl;
        }
    }
    
    cout << endl << "Shortage probabilities:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.blockingProbability[evaluatedQueue[i]] << endl;
    }
```

The model is finally evaluated by compiling the C++ program. If you have cloned the repository to your computer, remove the file `PythonWrapper.cpp`, and run `g++ -O3 *.cpp -o eval`. Run the program with `./eval`. 

### The complete example

```c++
#include "Model.h"

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

        //arrival rates for each customer type
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};
    
    //mean service time for each customer type
    vector<double> serviceTimes = {10,5,10,8};

    //capacity of each queue
    vector<int> capacity = {15,20,10,30};
    
    //fraction of rejected customers that are moved to an alternative queue node
    //this is an nCustomerTypes x nQueues matrix
    vector<vector<double>> relocationProbabilities = {{0.0,0.4,0.1,0.5},
                                                      {0.3,0.0,0.5,0.0},
                                                      {0.0,0.5,0.0,0.5},
                                                      {0.2,0.3,0.5,0.0}};
    
    //queue indices preferred by each customer type
    vector<int> preferredQueue = {0,1,2,3};

    //indices of the queues to be evaluated    
    vector<int> evaluatedQueue = {0,1,2,3};

   //create the model object
    Model mdl(arrivalRates,serviceTimes,capacity,
            relocationProbabilities,preferredQueue,
            evaluatedQueue);
    
    //now evaluate the model
    mdl.runModel();

    //get the results
    cout << "Occupancy distributions:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << "----" << "Queue " << evaluatedQueue[i] << "----" << endl;
        for (int j=0; j<mdl.queueDenDist[evaluatedQueue[i]].size(); j++){
            cout << mdl.queueDenDist[evaluatedQueue[i]][j] << endl;
        }
    }
    
    cout << endl << "Shortage probabilities:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.blockingProbability[evaluatedQueue[i]] << endl;
    }

    return 0;
}
```

## Model object and methods

### Model object

```c++
Model(vector<double> arrRates, //vector of arrival rates
        vector<double> serTimes, //vector of service times
        vector<int> cap, //vector of capacities
        vector<vector<double>> relProbMat, //relocation probability matrix
        vector<int> prefQ, //vector of preferred queues
        vector<int> evalQ, //queues to evaluate 
        string mdlt="auto", //sets the model type (auto, approximation or simulation)
        bool eqze=true); //If true, all service times are equalized and arrival rates adjusted to keep the load on the queues.
```

### Model methods

* `runModel()`. Evaluate the model.
* `setSeed(int sd)`. Set a seed.
* `setBurnIn(double bn)`. Set the burn-in time.
* `setMinimumSimulationTime(double mnTime)`. Set the simulation time. 
* `setMinSamples(int mnSamples)`. Set the minimum number of open/shortage samples.
* `setHyperStates(int openStates, int blockedStates)`. Set the number of phases in the hyper-exponential distributions accounting for the open/shortage time.


# Applications of RelSys

## Articles

Andersen, Anders Reenberg, Bo Friis Nielsen, and Andreas Lindhardt Plesner. 2023. "An Approximation of the Inpatient Distribution in Hospitals with Patient Relocation Using Markov Chains." Healthcare Analytics 3: 100145. https://doi.org/10.1016/j.health.2023.100145.

Leeters, Christoph, and Anders Reenberg Andersen. 2023. "Queueing systems with relocation of customers." ORbit: medlemsblad for Dansk Selskab for Operationsanalyse.

# How to cite

[![DOI](https://zenodo.org/badge/293829002.svg)](https://zenodo.org/badge/latestdoi/293829002)

Anders Reenberg Andersen. (2022). areenberg/RelSys: RelSys - first release (v1.0). Zenodo. https://doi.org/10.5281/zenodo.6037435

# License

Copyright 2023 Anders Reenberg Andersen.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
