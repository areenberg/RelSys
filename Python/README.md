# RelSys (**Rel**ocation **Sys**tem)
RelSys is a tool for evaluating a system of queues where arriving customers can be relocated to an alternative queue if none of the servers in the preferred queue are idle.

# Table of contents

1. Description of RelSys
2. How to use
3. How to cite
4. Licence

# Description of RelSys

Consider a number of parallel queues where the capacity of the queue equals the number of servers. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue (also denoted the Erlang loss or Erlang-B system), a customer is rejected and lost from the system if all the servers are occupied upon arrival. In RelSys, we allow customers to be transferred with a probability to one of the other queues. If there is an idle server in the alternative queue, the customer is accepted and served with an exponentially distributed time with the same rate-parameter as the customer would have in the preferred queue. The figure below depicts an example featuring two queues where customers are relocated (i.e. transferred) with a probability to the other queue whenever the preferred queue is full. 

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

## Installation

Download and install `relsys` directly from PyPI.

```
pip install relsys
```

## Usage

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

# Capsule on Code Ocean

We have published a capsule for the Linux CLI on Code Ocean. The URL and DOI for the capsule follows below:

* URL: https://codeocean.com/capsule/7104737/tree
* DOI: https://doi.org/10.24433/CO.2728562.v1

# How to cite

[![DOI](https://zenodo.org/badge/293829002.svg)](https://zenodo.org/badge/latestdoi/293829002)

# License

Copyright 2024 Anders Reenberg Andersen.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
