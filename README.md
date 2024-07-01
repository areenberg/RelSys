# RelSys (**Rel**ocation **Sys**tem)
RelSys is a tool for evaluating a system of queues where arriving customers can be relocated to an alternative queue if none of the servers in the preferred queue are idle.
We have developed two interfaces for RelSys: A Python module and a command-line interface.  

# Table of contents

1. Description of RelSys
   * Input parameters
   * Output types
2. How to use
   * Python
   * Command-line interface
   * C++
3. How to cite
4. Licence

# Description of RelSys

Consider a number of parallel queues where the capacity of the queue equals the number of servers. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue (also denoted the Erlang loss or Erlang-B system), a customer is rejected and lost from the system if all the servers are occupied upon arrival. In RelSys, we allow customers to be transferred with a probability to one of the other queues. If there is an idle server in the alternative queue, the customer is accepted and served with an exponentially distributed time with the same rate-parameter as the customer would have in the preferred queue. The figure below depicts an example featuring two queues where customers are relocated (i.e. transferred) with a probability to the other queue whenever the preferred queue is full. 

<img src="https://github.com/areenberg/RelSys/blob/development/images/example_system.jpeg" width="430" height="500">

## Input parameters

RelSys uses five types of input parameters:

* An arrival rate vector, where each element corresponds to a customer type.
* A service time vector, where each element corresponds to a customer type.
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
* Expected fraction of occupied capacity.

All outputs, except for the *expected occupancy* and *expected fraction of occupied capacity*, can be evaluated from two customer perspectives: The perspective of the customers preferring the queues, and the perspective of all customers arriving to the queues.    

# How to use

Below are guides on how to install and use the Python package and command-line interface.

## Python (3.10)

**Linux installation:**

Download and install `relsys` directly from PyPI.

```
pip install relsys
```

`relsys` works *out of the box* on Linux.


**Windows installation:**

Requires [Visual Studio 2022](https://visualstudio.microsoft.com/downloads/) (17.9) with MSVC C++ compiler and libraries installed (see below).

<img src="https://github.com/areenberg/RelSys/blob/master/images/visualstudioinstaller.png" alt="Logo" width="300">

After installing Visual Studio (incl. MSVC C++ compiler and libraries), install directly from PyPI with:

```
pip install relsys
```

**How to use the Python package:**

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

Create the model object and import the parameters,

```python
mdl = relsys.model(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)
```

The model can now be evaluated with `run`,

```python
mdl.run()
```

Return the resulting occupancy distributions with `getDensity` and shortage probabilities with `getShortageProb`,

```python
for queueIdx in range(4):
    print(mdl.getDensity(queueIdx))

for queueIdx in range(4):
    print(mdl.getShortageProb(queueIdx))
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

#create the model object and import the parameters
mdl = relsys.model(arrivalRates,serviceTimes,capacity,relocationProbabilities,preferredQueue)

#run the model
mdl.run()

#check the resulting occupancy distribution of each queue 
for queueIdx in range(4):
    print(mdl.getDensity(queueIdx))

#check the resulting shortage probabilities of each queue 
for queueIdx in range(4):
    print(mdl.getShortageProb(queueIdx))
```

## List of functions

### Model settings

* `setType(string mdltype)`. Set the method to use in the evaluation of the model ("simulation" (default), "approximation", "auto"). 
* `queuesEval(list qEvalIdx)`. Set the indices of queues to evaluate.
* `equalize(bool equalize)`. Specify if service times should be equalized and loads correspondingly adjusted (True: On, False: Off (default)).
* `setVerbose(bool set)`. Control verbose (True: On, False: Off (default)).
* `setSeed(int sd)`. Set the seed.
* `setAccSamType(string stype)`. Set the accuracy estimation type for the simulation ("preferred" (default), "all").
* `setSimTolerance(double tol)`. Set the tolerance level for the accuracy estimation in the simulation (default: 5e-3).
* `setBurnIn(double bin)`. Set the burn-in time of the simulation.
* `setSimTime(double mnTime)`. Set the simulation time.
* `setSamples(int mnSamples)`. Set the minimum number of open/shortage samples.
* `setHyperPhases(int openStates, int blockedStates)`. Set the number of phases in the hyper-exponential distributions accounting for the open/shortage time.

### Run calculations

* `run()`. Evaluate the model using the input parameters.

### Get results

* `getDensity(int queueIdx=0, string type="preferred")`. Return the density distribution of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getFreq(int queueIdx=0, string type="preferred")`. Return the frequency distribution of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getShortageProb(int queueIdx=0, string type="preferred")`. Return the shortage probability of a queue. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getAvailProb(int queueIdx=0, string type="preferred")`. Return the probability that at least one server is available. The second argument specifies the arrival type: "all" and "preferred" (default).
* `getExpOccupany(int queueIdx=0)`. Return the expected number of occupied servers.
* `getExpOccFraction(int queueIdx=0)`. Return the expected fraction of occupied servers.

### Return imported variables

* `getArrivalRates()`. Return the imported arrival rates.
* `getServiceTimes()`. Return the imported service times.
* `getCapacity()`. Return the imported capacities.
* `getReloc()`. Return the imported relocation probabilities.
* `getPreferredQueue()`. Return the imported preferred queues.

## Command-line Interface

We have created a Command-Line Interface (CLI), which is similar to the Python package in terms of features, inputs, and outputs. The CLI utilizes files to import the input parameters and export the results.

**Linux installation:**

Clone the repository:

```
git clone https://github.com/areenberg/RelSys.git
```

Navigate to the `RelSys` subfolder in the cloned repository and compile the EXE-file (`relsys.exe`) for the CLI: 

```
g++ -O3 -std=c++11 $(ls *.cpp | grep -v 'PythonWrapper.cpp\|ModuleInterface.cpp') -o relsys.exe
``` 

Make sure the CLI works: `./relsys.exe -demo`.

**Windows installation:**

Make sure that *g++* is installed and added to your *Path* variable:
1. To get *g++*: (1) Install [MinGW](https://sourceforge.net/projects/mingw/). (2) Start MinGW and install all basic support libraries for C++.
2. Add *C:\MinGW\bin* to the Windows *Path* variable.
3. Make sure *g++* works: `g++ --version`.

Now clone the repository:

```
git clone https://github.com/areenberg/RelSys.git
```

Navigate to the `RelSys` subfolder in the cloned repository, open *Windows PowerShell*, and compile the EXE-file (`relsys.exe`) for the CLI: 

```
g++ -O3 -std=c++11 -o relsys.exe Combinatorial.cpp LinSolver.cpp QueuePerformance.cpp Customer.cpp CustomerData.cpp main.cpp RelocEvaluation.cpp Model.cpp RelocSimulation.cpp EntireSystem.cpp Experiments.cpp PhaseFitter.cpp StatusBar.cpp HeuristicQueue.cpp SystemParameters.cpp HyperQueue.cpp QueueData.cpp
```

Make sure the CLI works: `./relsys.exe -demo`.

**How to use the CLI:**

The syntax for the CLI is `relsys [options]`. Use the `-help` flag to view all available options.

```
relsys -help
```

Create a space-separated file for each of the input parameter types. For instance,

*arrivalRates.txt*
```
0.8 2.5 0.6 2.8
```

*serviceTimes.txt*
```
10 5 10 8
```

*capacity.txt*
```
15 20 10 30
```

*relocProbs.txt*
```
0.0 0.4 0.1 0.5
0.3 0.0 0.5 0.0
0.0 0.5 0.0 0.5
0.2 0.3 0.5 0.0
```

*preferred.txt*
```
0 1 2 3
```

Evaluate the model using simulation, and save the result in a semicolon-separated file named `results.csv`.

```
relsys -m simulation -arr arrivalRates.txt -ser serviceTimes.txt -cap capacity.txt -rel relocProbs.txt -prq preferred.txt -o results.csv
```

## C++

The following guide will show you how to use the C++ source code for evaluating a model. Run `git clone https://github.com/areenberg/RelSys.git`. The directory `RelSys/` contains the complete source code for RelSys. Start by heading here. Replace the content of the current *main*-file with:

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
        for (int j=0; j<mdl.queueDenDistPref[evaluatedQueue[i]].size(); j++){
            cout << mdl.queueDenDistPref[evaluatedQueue[i]][j] << endl;
        }
    }
    
    cout << endl << "Shortage probabilities:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.blockingProbabilityPref[evaluatedQueue[i]] << endl;
    }
```

The model is finally evaluated by compiling the C++ program. If you have cloned the repository to your computer, remove the files `PythonWrapper.cpp` and `ModuleInterface.cpp`, and run `g++ -O3 -std=c++11 *.cpp -o relsys.exe` in your terminal. Run your program with `./relsys.exe`. 

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
        for (int j=0; j<mdl.queueDenDistPref[evaluatedQueue[i]].size(); j++){
            cout << mdl.queueDenDistPref[evaluatedQueue[i]][j] << endl;
        }
    }
    
    cout << endl << "Shortage probabilities:" << endl;
    for (int i=0; i<evaluatedQueue.size(); i++){
        cout << mdl.blockingProbabilityPref[evaluatedQueue[i]] << endl;
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
        string mdlt="simulation", //sets the model type (auto, approximation or simulation)
        bool eqze=false); //If true, all service times are equalized and arrival rates adjusted to keep the load on the queues.
```

### Model methods

* `runModel()`. Evaluate the model.
* `setSeed(int sd)`. Set a seed.
* `setBurnIn(double bn)`. Set the burn-in time.
* `setMinimumSimulationTime(double mnTime)`. Set the simulation time. 
* `setMinSamples(int mnSamples)`. Set the minimum number of open/shortage samples.
* `setHyperStates(int openStates, int blockedStates)`. Set the number of phases in the hyper-exponential distributions accounting for the open/shortage time.
* `setSimTolerance(double at)`. Set the tolerance for the automatic termination of the simulation.
* `setAccuracySampleType(string stype)`. Set the accuracy evaluation type for the automatic termination of the simulation (preferred , all).

# Applications of RelSys

## Articles

Andersen, Anders Reenberg, Bo Friis Nielsen, and Andreas Lindhardt Plesner. 2023. "An Approximation of the Inpatient Distribution in Hospitals with Patient Relocation Using Markov Chains." Healthcare Analytics 3: 100145. https://doi.org/10.1016/j.health.2023.100145.

Leeters, Christoph, and Anders Reenberg Andersen. 2023. "Queueing systems with relocation of customers." ORbit: medlemsblad for Dansk Selskab for Operationsanalyse.

# Capsule on Code Ocean

We have published a capsule for the Linux CLI on Code Ocean. The URL and DOI for the capsule follows below:

* URL: https://codeocean.com/capsule/7104737/tree
* DOI: https://doi.org/10.24433/CO.2728562.v1

# How to cite

Please cite our paper:

Andersen, Anders Reenberg, Bo Friis Nielsen, and Andreas Lindhardt Plesner. 2023. "An Approximation of the Inpatient Distribution in Hospitals with Patient Relocation Using Markov Chains." Healthcare Analytics 3: 100145. https://doi.org/10.1016/j.health.2023.100145.

or the software itself using the DOI below:

[![DOI](https://zenodo.org/badge/293829002.svg)](https://zenodo.org/badge/latestdoi/293829002)
