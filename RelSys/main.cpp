
/*
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 * Created on August 4, 2021, 02:11 PM
 */

#include "CustomerData.h"
#include "QueueData.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "SystemParameters.h"
#include "Model.h"

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;


int main(int argc, char** argv) {

    //--------------------------
    //PARSE ARGUMENTS
    //--------------------------

    if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0)
    {
        std::cout << " -f\tfile path of the model" << std::endl;
        std::cout << " -q\tqueue to simulate if using markov chain model. Use -1 to evaluate all queues." << std::endl;
        return 0;
    }

    char* filePath = NULL; // path of the model file
    int widx = -1; //queue index to be evaluated

    for(int i = 1; i < argc - 1; i++)
    {
        if(strcmp(argv[i], "-f") == 0)
        {
            filePath = argv[i+1];
            i++;
        }
        else if(strcmp(argv[i], "-q") == 0)
        {
            widx = atoi(argv[i+1]);
            i++;
        }
        else
        {
            std::cout << "Invalid argument: " << argv[i] << std::endl;
        }
    }

    if(filePath == NULL)
    {
        std::cout << "Missing argument: Please specify file path for the model using the -f flag" << std::endl;
    }

    //--------------------------
    //IMPORTING THE MODEL
    //--------------------------

    Model model;
    if(!model.ReadFromFile(filePath))
    {
        std::cout << "Error: Invalid model file path!" << std::endl;
        return 1;
    }

    if(strcmp(model.simulationMode, "markov") != 0 && strcmp(model.simulationMode, "simulation") != 0)
    {
        std::cout << "Error: Invalid simulation mode!" << std::endl;
        return 1;
    }

    //--------------------------
    //CUSTOMER SETUP
    //--------------------------

    //total number of customer types
    int nCustomerTypes = model.patientTypes.size();

    //create the customer type objects
    CustomerData * custs_array = new CustomerData[nCustomerTypes];
    for (int i=0; i<nCustomerTypes; i++){
        custs_array[i] = CustomerData(model.patientTypes[i].preferredWard,
                                      model.patientTypes[i].arrivalRate,
                                      model.patientTypes[i].meanLengthOfStay,
                                      model.patientTypes[i].relocationProbabilities);
    }

    //--------------------------
    //QUEUE SETUP
    //--------------------------

    //total number of queues
    int nQueues = model.wards.size();

    //calculate system input parameters from customer types to queues
    SystemParameters sysParam(nQueues,nCustomerTypes,custs_array);

    //now create the queue objects
    QueueData * wd_array = new QueueData[nQueues];
    for (int i=0; i<nQueues; i++){
        wd_array[i] = QueueData(i, //<<-- important: this input must always equal the index of the queue in the array
                                sysParam.queueArrivalRate(i),
                                sysParam.queueServiceRate(i),
                                model.wards[i].capacity,
                                sysParam.queueRelProbability(i));
    }

    if(strcmp(model.simulationMode, "markov") == 0)
    {
        //--------------------------
        //CTMC APPROX EVALUATION
        //--------------------------

        //setup model object
        RelocEvaluation mdl(nQueues,wd_array);

        //first simulate open/blocked time-windows
        int seed = 123;
        int bin = 365;
        int mt = 365;
        mdl.runSimulation(seed,bin,mt,50);
        
        if(widx == -1) // simulate all queues
        {
            for(int i = 0; i < nQueues; i++)
            {
                mdl.runHeuristic(i);

                if(i == 0)
                {
                    cout << "--- RESULTS ---" << endl;
                }

                cout << "--- queue " << i << " ---" << endl;

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
            }
        }
        else // simulate only one queue
        {
            mdl.runHeuristic(widx);
        }

        //get the result
        cout << "--- RESULTS ---" << endl;
        cout << "--- queue " << widx << " ---" << endl;

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
    }

    if(strcmp(model.simulationMode, "simulation") == 0)
    {
        //--------------------------
        //SIMULATION EVALUATION
        //--------------------------

        //setup simulation model object
        RelocSimulation sim_mdl(nQueues,wd_array);

        //setup and run simulation
        sim_mdl.setSeed(123); //set the seed
        double burnIn = 365; //burn-in time
        double minTime = 100000; //minimum simulation time
        vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
        sim_mdl.disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling

        sim_mdl.simulate(burnIn,minTime,maxWardSamples); //now run the simulation

        //In contrary to the heuristic evaluation, the entire system is evaluated at the
        //same time. However, for the sake of this demonstration we choose only to print
        //the results from a single queue.
        int sim_widx = 0;

        //get the result
        cout << "--- RESULTS ---" << endl;

        cout << "Marginal frequency distribution:" << endl;
        for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
            cout << sim_mdl.wardFreqDist[sim_widx][i] << endl;
        }
        cout << "Marginal density distribution:" << endl;
        for (int i=0; i<sim_mdl.wardFreqDist[sim_widx].size(); i++){
            cout << sim_mdl.wardDenDist[sim_widx][i] << endl; //corresponds to marginalDist in the heuristic evaluation
        }

        cout << "Probability of rejection:" << endl;
        cout << sim_mdl.blockingProbability[sim_widx] << endl;

        cout << "Expected server occupancy:" << endl;
        cout << sim_mdl.expectedOccupancy[sim_widx] << endl;

        cout << "Expected fraction of servers occupied:" << endl;
        cout << sim_mdl.expOccFraction[sim_widx] << endl;
    }

    return 0;
}
