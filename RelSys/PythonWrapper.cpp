#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>

#include "Model.h"

using namespace std;
namespace py = pybind11;


//TASKS TO GO
//- create on/off switch for verbose from model
//- create method for checking input parameters
//- move compiled library and example file to separate folder
//- commit and push to github and *update readme*

//input data structure
struct {
    vector<double> arrivalRates;
    vector<double> serviceTimes;
    vector<int> capacity;
    vector<vector<double>> relocationProbabilities;
    vector<int> preferredQueue;
} data;

//model settings structure
struct {
    bool evalAll=true;
    bool equalizeService=true;
    vector<int> evaluatedQueue;
    string modelType="auto";
} settings;

//results
struct {
    vector<double> shortageProbability;
    vector<double> availProbability;
    vector<vector<double>> queueDenDist;
    vector<vector<int>> queueFreqDist;
    vector<double> expectedOccupancy;
    vector<double> expOccFraction;
} results;


//import the data for the model
void importData(py::list arr, py::list ser, py::list cap,
py::list relProb, py::list prefQ){

    data.arrivalRates = arr.cast<vector<double>>();
    data.serviceTimes = ser.cast<vector<double>>();
    data.capacity = cap.cast<vector<int>>();
    data.relocationProbabilities = relProb.cast<vector<vector<double>>>();
    data.preferredQueue = prefQ.cast<vector<int>>();

}

bool parametersAvail(){
    if (!data.arrivalRates.empty()&&!data.serviceTimes.empty()&&
    !data.capacity.empty()&&!data.relocationProbabilities.empty()&&
    !data.preferredQueue.empty()){
        return(true);
    }else{
        return(false);    
    }
}

void evaluateAllQueues(){
    settings.evaluatedQueue.resize(data.capacity.size());
    for (int i=0; i<data.capacity.size(); i++){
        settings.evaluatedQueue[i] = i;
    }
}

void queuesToEvaluate(py::list qEvalIdx){
    //set the indices of queues to
    //evaluate
    settings.evalAll=false;
    settings.evaluatedQueue = qEvalIdx.cast<vector<int>>();
}

void setType(string mdltype){
    //set the type of model to use
    //in the evaluation
    settings.modelType = mdltype;
}

void equalizeService(bool equalize){
    //specifies if the service rates
    //in the model should be equalized
    //and loads correspondingly balanced
    settings.equalizeService = equalize;
}

py::list getArrivalRates(){
    return(py::cast(data.arrivalRates));
}

py::list getServiceTimes(){
    return(py::cast(data.serviceTimes));
}

py::list getCapacity(){
    return(py::cast(data.capacity));
}

py::list getRelocationProbabilities(){
    return(py::cast(data.relocationProbabilities));
}

py::list getPreferredQueue(){
    return(py::cast(data.preferredQueue));
}

py::list getDensityDistribution(int queueIndex){
    return(py::cast(results.queueDenDist[queueIndex]));
}

py::list getFrequencyDistribution(int queueIndex){
    return(py::cast(results.queueFreqDist[queueIndex]));
}

double getShortageProbability(int queueIndex){
    return(results.shortageProbability[queueIndex]);
}

double getAvailProbability(int queueIndex){
    return(results.availProbability[queueIndex]);
}

double getExpectedOccupancy(int queueIndex){
    return(results.expectedOccupancy[queueIndex]);
}

double getExpOccFraction(int queueIndex){
    return(results.expOccFraction[queueIndex]);
}

void runCalculations(){

    //--------------------------
    //MODEL DATA
    //--------------------------

    //check if input provided
    if(!parametersAvail()){
        cout << "The input parameters are not available. Please use importData() to import parameters for the model. Aborting program." << endl;
        exit(1);
    }
    //indices of queues to be evaluated by the model
    if (settings.evalAll){
        evaluateAllQueues();        
    }

    //check parameters are feasible
    //CODE HERE    
    
    //--------------------------
    //MODEL EVALUATION
    //--------------------------

    //create the model object
    Model mdl(data.arrivalRates,data.serviceTimes,data.capacity,
            data.relocationProbabilities,data.preferredQueue,
            settings.evaluatedQueue,settings.modelType,
            settings.equalizeService);
    
    //now evaluate the model
    mdl.runModel();
    
    //--------------------------
    // GET THE RESULTS
    //--------------------------

    //queue density distribution    
    results.queueDenDist.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueDenDist[settings.evaluatedQueue[i]].resize(mdl.queueDenDist[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueDenDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueDenDist[settings.evaluatedQueue[i]][j] = mdl.queueDenDist[settings.evaluatedQueue[i]][j];
        }
    }

    //queue frequency distribution    
    results.queueFreqDist.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueFreqDist[settings.evaluatedQueue[i]].resize(mdl.queueFreqDist[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueFreqDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueFreqDist[settings.evaluatedQueue[i]][j] = mdl.queueFreqDist[settings.evaluatedQueue[i]][j];
        }
    }
    
    //shortage probabilities
    results.shortageProbability.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.shortageProbability[settings.evaluatedQueue[i]] = mdl.blockingProbability[settings.evaluatedQueue[i]];
    }

    //availability probabilities
    results.availProbability.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.availProbability[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbability[settings.evaluatedQueue[i]];
    }

    //expected occupancy
    results.expectedOccupancy.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expectedOccupancy[settings.evaluatedQueue[i]] = mdl.expectedOccupancy[settings.evaluatedQueue[i]];
    }

    //expected fraction of customers occupied
    results.expOccFraction.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expOccFraction[settings.evaluatedQueue[i]] = mdl.expOccFraction[settings.evaluatedQueue[i]];
    }
}

PYBIND11_MODULE(relsys, m) {
    
    //import data
    m.def("input",&importData,"Import the input data to the model.");

    //model settings
    m.def("setType",&setType,"Specify the method to use in the evaluation of the model (auto, simulation, approximation)."); 
    m.def("queuesEval",&queuesToEvaluate,"Specify indices of queues to evaluate.");
    m.def("equalize",equalizeService,"Specify if service times should be equalized and loads correspondingly adjusted (True=On, False=Off).");

    //run calculations    
    m.def("run", &runCalculations,"Evaluate the model using the input parameters.");
    
    //return results
    m.def("getDensity",&getDensityDistribution,"Return the density distribution of a queue.");
    m.def("getFreq",&getFrequencyDistribution,"Return the frequency distribution of a queue.");
    m.def("getShortageProb",&getShortageProbability,"Return the shortage probability of a queue.");
    m.def("getAvailProb",&getAvailProbability,"Return the probability that at least one server is available.");
    m.def("getExpOccupany",&getExpectedOccupancy,"Return the expected number of occupied servers.");
    m.def("getExpOccFraction",&getExpOccFraction,"Return the expected fraction of occupied servers.");

    //return imported variables  
    m.def("getArrivalRates",&getArrivalRates,"Return the imported arrival rates.");
    m.def("getServiceTimes",&getServiceTimes,"Return the imported service times.");
    m.def("getCapacity",&getCapacity,"Return the imported capacities.");
    m.def("getReloc",&getRelocationProbabilities,"Return the imported relocation probabilities.");
    m.def("getPreferredQueue",&getPreferredQueue,"Return the imported preferred queues.");
        
}
