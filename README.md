# RelSys
RelSys is a tool for evaluating a system of queues with finite capacity and multiple classes of customers. The queues are connected, but *only* through customer relocations (i.e. customers that are transferred to an alternative queue if they are blocked instead of being rejected from the system).
The tool is written in C++ and currently employs two different approaches for evaluating the system. The first (`RelocEvaluation`) is a mathematical model and is based on a continuous-time Markov chain approximation of the system, and the second (`RelocSimulation`) is a discrete-event simulation. The input is the same in both cases, but each approach comes with different methods for evaluating the system and viewing the results.

# Table of contents

1. How does it work
2. Getting started
3. How to cite
4. Licence

# How does it work

Consider a number of parallel queues with finite capacity. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue (also denoted the Erlang loss or Erlang-B system), a customer is rejected and lost from the system if all the servers are occupied upon arrival. However, in the RelSys-modeling tool, we allow customers to be transferred with a probability to one of the other queues. If there is an idle server in the alternative queue, the customer is accepted and served with an exponentially distributed time with the same rate-parameter as the customer would have had in the original queue. The figure below depicts an example featuring two queues where customers are relocated (i.e. transferred) with a probability to the other queue whenever the preferred queue is full. 

<img src="https://github.com/areenberg/RelSys/blob/development/images/example_system.jpeg" width="430" height="500">

## Briefly about the setup

The RelSys tool is divided into:

1. *Customer setup*, where objects for each customer type are created using the `CustomerData` class.
2. *Queue setup*, where objects for each queue are created using the `QueueData` class. Data from the customer objects are used in the creation of the queue objects through a class called `SystemParameters` (see the figure below).
3. *Model choice and evaluation*, where a modeling approach is selected (either `RelocEvaluation` or `RelocSimulation`) and a solution is found.
4. *Get the results*, where results are printed/saved using the methods from the model objects. E.g. `mdl.blockingProbability` to get the probability of customer blocking.

The figure below shows how the various classes and input parameters are connected. 

<img src="https://github.com/areenberg/RelSys/blob/master/images/programSetupGraph.jpg" width="500" height="650">

## Input parameters

There are seven types of input parameters:

* Arrival rates (`double`). Used in the `CustomerData` class.
* Service times (`double`). Used in the `CustomerData` class.
* Relocation probabilities (`vector<double>`). Used in the `CustomerData` class.
* The queue preferred by the customer type (`int`). Used in the `CustomerData` class.
* Capacities (`int`). Used in the `QueueData` class.
* Number of customer types (`int`). Used in the `SystemParameters` class, and to create an array of `CustomerData` objects.
* Number of queues (`int`). Used in `SystemParameters`, the model class `RelocEvaluation`/`RelocSimulation`, and to create an array of `QueueData` objects.  

The relation between these parameters and the various classes is depicted with grey arrows in the figure above.

## Run and get results (Output methods)

As previously mentioned, the system can be evaluated using one of two modeling approaches:

1. Using a continuous-time Markov chain (CTMC) approximation of the system through `RelocEvaluation`.
2. Using a discrete-event simulation through `RelocSimulation`.

Both approaches take the number of queues and an array of `QueueData` objects, e.g. `RelocEvaluation mdl(nQueues,wd_array)`. The subsequent setup varies depending on the class.

## `RelocEvaluation`

After creating the model object, and prior to evaluating the system, a simulation of the internal dynamics have to be conducted. This is done using the `runSimulation(int seed, int burnIn, int minTime, int minSamples)` method. Here, `int seed` is the simulation seed, `int burnIn` is the burn-in time, `int minTime` is the minimum simulation-time the simulation has to run (this can often be set equal to `burnIn`) and `int minSamples` is the minimum number of samples that are collected by the simulation. This parameter has a default value of 50.

Now, the system can be evaluated using the method `runHeuristic(int widx)`, where `int widx` is the index of the queue being evaluated.

### Getting the results

The `RelocEvaluation` class features four types of performance measures. These are:

* Marginal probability distribution: `vector<double> marginalDist`.
* Probability of rejection: `double blockingProbability`.
* Expected server occupancy: `double expectedOccupancy`.
* Expected fraction of servers occupied: `double expOccFraction`. 

## `RelocSimulation`

Similar to `RelocEvaluation` we need to conduct some preliminary setup before the system can be evaluated. First, a seed is chosen using the `setSeed(int seed)` method. Next, we disable some of the features that we dont need using the `disableTimeSampling()` method. This speeds up the simulation.

Now, the system can be evaluated using `simulate(double burnIn, double minTime, vector<int> maxWardSamples, int minSamples)`, where `double burnIn` is the burn-in time, `double minTime` is the minimum simulation time, `vector<int> maxWardSamples` is an upper bound on the number of samples the simulation is allowed to collect from each queue. In most cases, this functionality should be disabled, which is done by setting the input vector to `vector<int> maxWardSamples(1,-1)` (i.e. with one element that has a value of -1). The final parameter, `int minSamples`, is the minimum number of samples of periods where the queues are either open or in shortage of capacity. This parameter has a default value of 50, but has no effect on the results if the method `disableTimeSampling()` has been used.

### Getting the results

The `RelocSimulation` class features five types of performance measures. These are:

* Marginal frequency distribution: `vector<vector<int>> wardFreqDist`.
* Marginal probability distribution: `vector<vector<double>> wardDenDist`.
* Probability of rejection: `vector<double> blockingProbability`.
* Expected server occupancy: `vector<double> expectedOccupancy`.
* Expected fraction of servers occupied: `vector<double> expOccFraction`. 

# Getting started

Start by loading the fundamental classes.

```c++

#include "CustomerData.h" //information about each customer type
#include "SystemParameters.h" //information about the overall system parameters
#include "QueueData.h" //information about each queue
#include "RelocEvaluation.h" //to evaluate the system using the CTMC approximation
#include "RelocSimulation.h" //to evaluate the system using discrete-event simulation

```

Now, specify the customer type parameters. In this example we consider a system of four types.

```c++

    //total number of customer types
    int nCustomerTypes = 4;

    //arrival rates for each customer type
    vector<double> arrivalRates = {1.2,1.8,0.5,2.5};

    //mean service time for each customer type
    vector<double> serviceTimes = {10,50,25,5};

    //fraction of rejected customers that are moved to an alternative queue node
    //this is an nCustomerTypes x nQueues matrix
    vector<vector<double>> relProbs = {{0.0,0.4,0.1,0.5},
                                       {0.0,0.0,0.2,0.3}, //<<-- note: these do not have to sum to one
                                       {0.2,0.0,0.0,0.8},
                                       {0.9,0.05,0.05,0.0}};

    //queue indices preferred by each customer type
    vector<int> preferredQueue = {0,1,2,3};

    //create the customer type objects
    CustomerData * custs_array = new CustomerData[nCustomerTypes];
    for (int i=0; i<nCustomerTypes; i++){
        custs_array[i] = CustomerData(preferredQueue[i],
                                      arrivalRates[i],
                                      serviceTimes[i],
                                      relProbs[i]);
    }

```

Note that the relocation probabilities (`relProbs`) do not have to sum to one. If the sum is less than one a fraction of the customers are lost without trying to relocate them first. The next step is to create the queues.
    
```c++    
    
    //total number of queues
    int nQueues = 4;

    //capacity of each queue
    vector<int> capacity = {15,18,12,16};

    //calculate system input parameters from customer types to queues
    SystemParameters sysParam(nQueues,nCustomerTypes,custs_array);

    //now create the queue objects
    QueueData * wd_array = new QueueData[nQueues];
    for (int i=0; i<nQueues; i++){
        wd_array[i] = QueueData(i,  //<<-- important: this parameter must not deviate from the order of the queues in the array
                                sysParam.queueArrivalRate(i),
                                sysParam.queueServiceRate(i),
                                capacity[i],
                                sysParam.queueRelProbability(i));
    }

```

The queueing objects are now ready to be plugged into the model object. In the following, we show how to create a model object, evaluate the system and get the results using both the heuristic approach and the discrete-event simulation, respectively. You will notice that the procedure is almost the same in both cases.

## `RelocEvaluation` - Evaluation with a CTMC approximation

First, create the model object using the aforementioned queueing array, `wd_array`.

```c++

    //setup model object
    RelocEvaluation mdl(nQueues,wd_array);

```

Next, some minor setup is required. Important: The CTMC approach evaluates only *one queue at a time*. Thus, the target queue needs to be specified. In the following, we do this using the variable `widx`.

```c++

    //minor setup and preliminary calculations
    int seed = 123;
    int bin = 365;
    int minTime = 365;
    int minSamples = 50
    mdl.runSimulation(seed,bin,minTime,minSamples);
    
    //choose a queue to evaluate
    int widx = 0; //queue index to be evaluated

```

The system is now ready to be evaluated.

```c++

    mdl.runHeuristic(widx);

```

Note that the memory consumption as well as the runtime might be *very* high depending on the size of the system.

All four performance measures are stored in the model object when the evaluation of the system is complete. 

```c++

    cout << "--- RESULTS ---" << endl;
   
    cout << "Marginal probability distribution:" << endl;
    for (int i=0; i<mdl.marginalDist.size(); i++){
        cout << mdl.marginalDist[i] << endl;
    }
    
    cout << "Probability of rejection:" << endl;
    cout << mdl.blockingProbability << endl;

    cout << "Expected server occupancy:" << endl;
    cout << mdl.expectedOccupancy << endl;

    cout << "Expected fraction of servers occupied:" << endl;
    cout << mdl.expOccFraction << endl;

```

## `RelocSimulation` - Evaluation by simulation

Once again, create the model object using the aforementioned queueing array, `wd_array`.

```c++

    //setup simulation model object
    RelocSimulation sim_mdl(nQueues,wd_array);
    
```

And again, some minor setup is required before we can proceed. The simulation evaluates all queues at the same time, so this time we do not have to settle on a single queue.

```c++

    //setup and run simulation
    sim_mdl.setSeed(123); //set the seed
    double burnIn = 365; //burn-in time
    double minTime = 50000; //minimum simulation time
    vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
    sim_mdl.disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling
    
```

The system can now be evaluated.

```c++

    sim_mdl.simulate(burnIn,minTime,maxWardSamples);

```

The simulation has much lower memory requirements than the CTMC approach and therefore able to evaluate much larger systems. The runtime and the precision depends greatly on when the simulation is terminated. In the above, the minimum simulation-time is specified using `minTime`. 

Recall that in contrary to the heuristic approach, the simulation evaluates the entire system at the same time. However, for the sake of this demonstration we chose only to print the results from a single queue, `sim_widx`.

```c++
    
    int sim_widx = 0;

    cout << "Marginal frequency distribution:" << endl;
    for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
        cout << sim_mdl.wardFreqDist[sim_widx][i] << endl;
    }
    cout << "Marginal probability distribution:" << endl;
    for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
        cout << sim_mdl.wardDenDist[sim_widx][i] << endl;
    }
    
    cout << "Probability of rejection:" << endl;
    cout << sim_mdl.blockingProbability[sim_widx] << endl;

    cout << "Expected server occupancy:" << endl;
    cout << sim_mdl.expectedOccupancy[sim_widx] << endl;

    cout << "Expected fraction of servers occupied:" << endl;
    cout << sim_mdl.expOccFraction[sim_widx] << endl;

```

# How to cite

Coming soon.

# License

Copyright 2021 Anders Reenberg Andersen.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
