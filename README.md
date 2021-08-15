# RelSys (Not up to date)
RelSys is a tool for evaluating a system of queues with finite capacity and multiple classes of customers. The queues are connected, but *only* through customer relocations (i.e. customers that are transferred to an alternative queue instead of being rejected from the system).
The tool is written in C++ and currently employs two different approaches for evaluating the system. The first (`RelocEvaluation`) is a heuristic mathematical model, which is based on a continuous-time Markov chain, and the second (`RelocSimulation`) is a discrete-event simulation. The input is the same in both cases, but each approach comes with different methods for evaluating the system and viewing the results.

# Table of contents

1. How does it work
2. Getting started
3. Further details
4. How to cite
5. Licence

# How does it work

Consider a number of parallel queues with finite capacity. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue (also denoted the Erlang loss or Erlang-B system), a customer is rejected and lost from the system if all the servers are occupied upon arrival. However, in the RelSys-modeling tool, we allow customers to be transferred with a probability to one of the other queues. If there is an idle server in the alternative queue, the customer is accepted and served with an exponentially distributed time with the same rate-parameter as the customer would have had in the original queue. Thus the system is multiclass. The figure below depicts an example featuring two queues where customers are relocated (i.e. transferred) with a probability to the other queue whenever the preferred queue is full. 

<img src="https://github.com/areenberg/RelSys/blob/development/images/example_system.jpeg" width="400" height="500">

## Briefly about the setup

The setup is divided into:

1. *Customer setup*, where objects for each customer type are created using the `CustomerData` class.
2. *Queue setup*, where objects for each queue are created using the `QueueData` class. Data from the customer objects are used in the creation of the queue objects through a class called `SystemParameters`.
3. *Model choice and evaluation*, where a modeling approach is selected (either `RelocEvaluation` or `RelocSimulation`) and a solution is found.
4. *Get the results*, where results are printed/saved using the methods from the model objects. E.g. `mdl.blockingProbability` to get the probability of customer blocking.


<img src="https://github.com/areenberg/RelSys/blob/master/images/programSetupGraph.jpg" width="550" height="600">

## Input parameters

There are four overall input parameters:

* Arrival rates.
* Service rates.
* Capacities.
* Relocation probabilities (a vector for each queue).

# Getting started

Start by loading the fundamental classes.

```c++

#include "QueueData.h" //information about each queue
#include "RelocEvaluation.h" //to evaluate the system using the heuristic approach
#include "RelocSimulation.h" //to evaluate the system using discrete-event simulation

```

Now, specify the parameters. In this example we consider a system of four queues.

```c++

    //total number of queues    
    int nQueues = 4; 
    
    //arrival rates for each queue
    vector<double> arrivalRates = {1.0,0.5,2.5,2.0};
    
    //service rates for each queue
    vector<double> serviceRates = {0.04,0.02,0.1,0.1};
    
    //capacity of each queue
    vector<int> capacity = {30,40,10,15};
    
    //fraction of rejected customers that are moved to an alternative queue node
    vector<vector<double>> relProbs = {{0.0,0.1,0,8,0.1},
                                       {0.1,0.0,0.3,0.1},
                                       {0.4,0.5,0.0,0.1},
                                       {0.2,0.5,0.3,0.0}};

```

Note that the relocation probabilities do not have to sum to one. If the sum is less than one a fraction of the customers are lost without trying to relocate them first. The next step is to create the queues.
    
```c++    
    
    //the queueing objects
    QueueData * wd_array = new QueueData[nQueues];
    wd_array[0] = QueueData(0,arrivalRates[0],serviceRates[0],capacity[0],relProbs[0]);
    wd_array[1] = QueueData(1,arrivalRates[1],serviceRates[1],capacity[1],relProbs[1]);
    wd_array[2] = QueueData(2,arrivalRates[2],serviceRates[2],capacity[2],relProbs[2]);
    wd_array[3] = QueueData(3,arrivalRates[3],serviceRates[3],capacity[3],relProbs[3]);
    
```

The queueing objects are now ready to be plugged into the model object. In the following, we show how to create a model object, evaluate the system and get the results using both the heuristic approach and the discrete-event simulation, respectively. You will notice that the procedure is almost the same in both cases.

## `RelocEvaluation` - Heuristic evaluation

First, create the model object using the aforementioned queueing array, `wd_array`.

```c++

    //setup model object
    RelocEvaluation mdl(nQueues,wd_array);

```

Next, some minor setup is required. Important: The heuristic approach evaluates only *one queue at a time*. Thus, the target queue needs to be specified. In the following, we do this using the variable `widx`.

```c++

    //minor setup and preliminary calculations
    int seed = 123;
    double bin = 365;
    double minTime = 365;
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

### Getting the results

The `RelocEvaluation` class features four types of performance measures. These are:

* Marginal probability distribution: `vector<double> marginalDist`.
* Probability of rejection: `double blockingProbability`.
* Expected server occupancy: `double expectedOccupancy`.
* Expected fraction of servers occupied: `double expOccFraction`. 

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

The simulation has much lower memory requirements than the heuristic approach and therefore able to evaluate much larger systems. The runtime and the precision depends greatly on when the simulation is terminated. In the above, the minimum simulation-time is specified using `minTime`. 

### Getting the results

The `RelocSimulation` class features five types of performance measures. These are:

* Marginal frequency distribution: `vector<vector<int>> wardFreqDist`.
* Marginal probability distribution: `vector<vector<double>> wardDenDist`.
* Probability of rejection: `vector<double> blockingProbability`.
* Expected server occupancy: `vector<double> expectedOccupancy`.
* Expected fraction of servers occupied: `vector<double> expOccFraction`. 

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

# Further details

Further details about the initial setup and methods in `RelocEvaluation` and `RelocSimulation`.
Coming soon. 

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
