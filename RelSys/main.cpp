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

#include <Windows.h> // temp

using namespace std;

bool use_cli_visuals = true;

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
    bool forceWardIndex = false;

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
            forceWardIndex = true;
            i++;
        }
        else if(strcmp(argv[i], "-hide") == 0)
        {
            use_cli_visuals = false;
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
    cout << "#Importing the model" << endl;

    Model model;
    if(!model.ReadFromFile(filePath))
    {
        std::cout << "Error: Invalid model file path!" << std::endl;
        return 1;
    }

    cout << "#Model imported successfully" << endl;

    if(strcmp(model.simulationMode, "markov") != 0 && strcmp(model.simulationMode, "simulation") != 0)
    {
        std::cout << "Error: Invalid simulation mode!" << std::endl;
        return 1;
    }

    if(!forceWardIndex)
    {
        widx = model.wardIndex;
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
        cout << "customer[" << i << "] pw: " << model.patientTypes[i].preferredWard << ", ar: " << model.patientTypes[i].arrivalRate << ", ls: " << model.patientTypes[i].meanLengthOfStay << ", rp: {";
        for(int j = 0; j < model.patientTypes[i].relocationProbabilities.size(); j++){
            cout << model.patientTypes[i].relocationProbabilities[j] << ", ";
        }
        cout << "}" << endl;
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
        int seed = model.seed;
        int bin = model.burnInTime;
        int mt = model.minTime;
        int ms = model.minSamples;
        int ohs = model.openHyperStates;
        int bhs = model.blockedHyperStates;

        cout << "#seed: " << model.seed << endl;
        cout << "#burn-in: " << model.burnInTime << endl;
        cout << "#min time" << model.minTime << endl;
        cout << "#min samples" << model.minSamples << endl;
        cout << "#open blocks" << model.openHyperStates << endl;
        cout << "#blocked blocks" << model.blockedHyperStates << endl;

        cout << "#Starting in approximation mode" << endl;
        mdl.runSimulation(seed,bin,mt,ms);

        //choose number of states in PH distributions (optional)
        mdl.setOpenHyperStates(ohs);
        mdl.setBlockedHyperStates(bhs);
        
        if(widx == -1) // simulate all queues
        {
            cout << "--- RESULTS ---" << endl;
            for(int i = 0; i < nQueues; i++)
            {
                mdl.runHeuristic(i);

                cout << "--- queue " << i << " ---" << endl;

                cout << "Marginal density distribution:" << endl;
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

            //get the result
            cout << "--- RESULTS ---" << endl;
            cout << "--- queue " << widx << " ---" << endl;

            cout << "Marginal density distribution:" << endl;
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
    else
    {
        //--------------------------
        //SIMULATION EVALUATION
        //--------------------------

        //setup simulation model object
        RelocSimulation sim_mdl(nQueues,wd_array);

        //setup and run simulation
        sim_mdl.setSeed(model.seed); //set the seed
        double burnIn = model.burnInTime; //burn-in time
        double minTime = model.minTime; //minimum simulation time
        vector<int> maxWardSamples(1,-1); //disables the limit on occupancy samples
        sim_mdl.disableTimeSampling(); //speed-up the simulation by disabling the open/blocked time-window sampling

        cout << "#seed: " << model.seed << endl;
        cout << "#burn-in: " << model.burnInTime << endl;
        cout << "#min time" << model.minTime << endl;

        cout << "#Starting in simulation mode" << endl;
        sim_mdl.simulate(burnIn,minTime,maxWardSamples); //now run the simulation

        //get the result
        cout << "--- RESULTS ---" << endl;

        for (int sim_widx = 0; sim_widx < nQueues; sim_widx++){
            cout << "--- queue " << sim_widx << " ---" << endl;
            
            cout << "Marginal frequency distribution:" << endl;
            for(int j = 0; j < sim_mdl.wardFreqDist[sim_widx].size(); j++){
                cout << sim_mdl.wardFreqDist[sim_widx][j] << endl;
            }
        

            cout << "Marginal density distribution:" << endl;
            for(int j = 0; j < sim_mdl.wardDenDist[sim_widx].size(); j++){
                cout << sim_mdl.wardDenDist[sim_widx][j] << endl; //corresponds to marginalDist in the heuristic evaluation
            }
        

            cout << "Probability of rejection:" << endl;
            cout << sim_mdl.blockingProbability[sim_widx] << endl;

            cout << "Expected server occupancy:" << endl;
            cout << sim_mdl.expectedOccupancy[sim_widx] << endl;

            cout << "Expected fraction of servers occupied:" << endl;
            cout << sim_mdl.expOccFraction[sim_widx] << endl;
        }
    }

    cout << "--- END ---";
    return 0;
}