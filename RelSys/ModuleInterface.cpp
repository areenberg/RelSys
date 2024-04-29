
#include "ModuleInterface.h"

ModuleInterface::ModuleInterface() {
}

ModuleInterface::ModuleInterface(const ModuleInterface& orig) {
}

ModuleInterface::~ModuleInterface() {
}

//import the data for the model
void ModuleInterface::importData(py::list arr, py::list ser, py::list cap,
py::list relProb, py::list prefQ){

    data.arrivalRates = arr.cast<vector<double>>();
    data.serviceTimes = ser.cast<vector<double>>();
    data.capacity = cap.cast<vector<int>>();
    data.relocationProbabilities = relProb.cast<vector<vector<double>>>();
    data.preferredQueue = prefQ.cast<vector<int>>();

}

void ModuleInterface::setAccuracySampleType(string stype){
    settings.accSampleType = stype;    
}

void ModuleInterface::setSimulationTolerance(double tol){
    settings.accTol = tol;
}

void ModuleInterface::setSeed(int sd){
    settings.seed = sd;
}

void ModuleInterface::setBurnIn(double bn){
    if (bn>0){
        settings.burnIn=bn;
    }else{
        py::print("Burn-in time must be larger than 0. Aborting program.");
        exit(1);
    }
}

void ModuleInterface::setMinimumSimulationTime(double mnTime){
    if (mnTime>0){
        settings.minimumSimulationTime=mnTime;
    }else{
        py::print("Simulation time must be larger than 0. Aborting program.");
        exit(1);
    }
}

void ModuleInterface::setMinSamples(int mnSamples){
    if (mnSamples>0){
        settings.minSamples=mnSamples;
    }else{
        py::print("Min. number of open/shortage samples must be larger than 0. Aborting program.");
        exit(1);
    }
}

void ModuleInterface::setHyperStates(int openStates, int blockedStates){
    if (openStates>0 && blockedStates>0){
        settings.hyperStates[0] = openStates;
        settings.hyperStates[1] = blockedStates;
    }else{
        py::print("The number of open and blocked states must be larger than 0. Aborting program.");
        exit(1);
    }
}

bool ModuleInterface::parametersAvail(){
    if (!data.arrivalRates.empty()&&!data.serviceTimes.empty()&&
    !data.capacity.empty()&&!data.relocationProbabilities.empty()&&
    !data.preferredQueue.empty()){
        return(true);
    }else{
        return(false);    
    }
}

void ModuleInterface::setVerbose(bool set){
    settings.verbose=set;
}

void ModuleInterface::checkParameters(){
    //checks if parameters are feasible

    if (data.arrivalRates.size() != data.serviceTimes.size()){
        py::print("The number of arrival rates must equal the number of service times. Aborting program.");
        exit(1);
    }
    if (data.arrivalRates.size() != data.relocationProbabilities.size()){
        py::print("The number of arrival rates, and service times, must equal the number of rows in the relocation matrix. Aborting program.");
        exit(1);
    }
    if (data.arrivalRates.size() != data.preferredQueue.size()){
        py::print("The number of arrival rates, and service times, must equal the length of the vector specifying the preferred queues. Aborting program.");
        exit(1);
    }
    if (settings.evaluatedQueue.size()>data.capacity.size() || settings.evaluatedQueue.empty()){
        py::print("The length of the vector of evaluated queues must be between 1 and n_queues. Aborting program.");
        exit(1);
    }        
    for (int i=0; i<data.relocationProbabilities.size(); i++){
        if (data.capacity.size() != data.relocationProbabilities[i].size()){
            py::print("The length of the capacity vector must equal the number of columns in the relocation matrix. Aborting program.");
            exit(1);
        }
    }        
    for (int i=0; i<data.arrivalRates.size(); i++){
        if (data.arrivalRates[i]<=0.0 || data.serviceTimes[i]<=0.0 ){
            py::print("Arrival rates and service times must be larger than 0.0. Aborting program.");
            exit(1);
        }
        for (int j=0; j<data.relocationProbabilities[i][j]; j++){
            if (data.relocationProbabilities[i][j]<0.0 || data.relocationProbabilities[i][j]>1.0){
                py::print("Values in the relocation matrix must be between 0.0 and 1.0. Aborting program.");
                exit(1);
            }
        }
    }
    for (int i=0; i<data.preferredQueue.size(); i++){
        if (data.preferredQueue[i]<0 || data.preferredQueue[i]>(data.capacity.size()-1)){
            py::print("Indices in the vector of preferred queues must be between 0 and n_queues-1. Aborting program.");
            exit(1);
        }
    }
    for (int i=0; i<data.capacity.size(); i++){
        if (data.capacity[i]<1){
            py::print("The capacity of each queue must be equal to or larger than 1. Aborting program.");
            exit(1);
        }
    }          
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        if (settings.evaluatedQueue[i]<0 || settings.evaluatedQueue[i]>(data.capacity.size()-1)){
            py::print("Indices in the vector of evaluated queues must lie in the interval between 0 and n_queues-1. Aborting program.");
            exit(1);
        }    
    }
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        for (int j=0; j<settings.evaluatedQueue.size(); j++){
            if (i!=j && settings.evaluatedQueue[i]==settings.evaluatedQueue[j]){
                py::print("All indices in the vector of evaluated queues must be unique. The same queue cannot be evaluated twice. Aborting program.");
                exit(1);
            }
        }    
    }
    double sm;
    for (int i=0; i<data.relocationProbabilities.size(); i++){
        sm=0;
        for (int j=0; j<data.relocationProbabilities[i].size(); j++){
            sm+=data.relocationProbabilities[i][j];
        }
        if (sm>1.0){
            string out = "The sum of the relocation probabilities in row "+to_string(i)+" is equal to "+to_string(sm)+". The sum must be equal to or smaller than 1.0. Aborting program.";
            py::print(out);
            exit(1);
        }
    }
    if ((settings.minimumSimulationTime!=-1 && settings.burnIn==-1) || (settings.minimumSimulationTime==-1 && settings.burnIn!=-1)){
        py::print("It is not possible to set the overall simulation time without setting the burn-in time and vice versa. Aborting program.");
        exit(1);
    }
    if (settings.minimumSimulationTime!=-1 && settings.burnIn!=-1 && settings.minimumSimulationTime<=settings.burnIn){
        py::print("The simulation time has to be longer than the burn-in time. Aborting program.");
        exit(1);
    }
    if (settings.accTol<=0.0 || settings.accTol>=1.0){
        py::print("The tolerance level for the accuracy estimation procedure must be between 0 and 1 (e.g. 1e-2). Aborting program.");
        exit(1);
    }        

}

void ModuleInterface::evaluateAllQueues(){
    settings.evaluatedQueue.resize(data.capacity.size());
    for (int i=0; i<data.capacity.size(); i++){
        settings.evaluatedQueue[i] = i;
    }
}

void ModuleInterface::queuesToEvaluate(py::list qEvalIdx){
    //set the indices of queues to
    //evaluate
    settings.evalAll=false;
    settings.evaluatedQueue = qEvalIdx.cast<vector<int>>();
}

void ModuleInterface::setType(string mdltype){
    //set the type of model to use
    //in the evaluation
    settings.modelType = mdltype;
}

void ModuleInterface::equalizeService(bool equalize){
    //specifies if the service rates
    //in the model should be equalized
    //and loads correspondingly balanced
    settings.equalizeService = equalize;
}

py::list ModuleInterface::getArrivalRates(){
    if (parametersAvail()){
        return(py::cast(data.arrivalRates));
    }else{
        py::print("The requested parameter has not been imported. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getServiceTimes(){
    if (parametersAvail()){
        return(py::cast(data.serviceTimes));
    }else{
        py::print("The requested parameter has not been imported. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getCapacity(){
    if (parametersAvail()){
        return(py::cast(data.capacity));
    }else{
        py::print("The requested parameter has not been imported. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getRelocationProbabilities(){
    if (parametersAvail()){
        return(py::cast(data.relocationProbabilities));
    }else{
        py::print("The requested parameter has not been imported. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getPreferredQueue(){
    if (parametersAvail()){
        return(py::cast(data.preferredQueue));
    }else{
        py::print("The requested parameter has not been imported. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getDensityDistribution(int queueIndex, string type){
    if (results.evaluated){
        if (type.compare("all")==0){
            return(py::cast(results.queueDenDist[queueIndex]));
        }else if (type.compare("preferred")==0){
            return(py::cast(results.queueDenDistPref[queueIndex]));
        }else{
            py::print("Output type not recognized. Choose between 'all' and 'preferrred'. Aborting program.");
            exit(1);    
        }
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

py::list ModuleInterface::getFrequencyDistribution(int queueIndex, string type){
    if (results.evaluated){
        if (type.compare("all")==0){
            return(py::cast(results.queueFreqDist[queueIndex]));
        }else if (type.compare("preferred")==0){
            return(py::cast(results.queueFreqDistPref[queueIndex]));
        }else{
            py::print("Output type not recognized. Choose between 'all' and 'preferrred'. Aborting program.");
            exit(1);    
        }
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

double ModuleInterface::getShortageProbability(int queueIndex, string type){
    if (results.evaluated){
        if (type.compare("all")==0){
            return(results.shortageProbability[queueIndex]);
        }else if (type.compare("preferred")==0){
            return(results.shortageProbabilityPref[queueIndex]);
        }else{
            py::print("Output type not recognized. Choose between 'all' and 'preferrred'. Aborting program.");
            exit(1);    
        }
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

double ModuleInterface::getAvailProbability(int queueIndex, string type){
    if (results.evaluated){
        if (type.compare("all")==0){
            return(results.availProbability[queueIndex]);
        }else if (type.compare("preferred")==0){
            return(results.availProbabilityPref[queueIndex]);
        }else{
            py::print("Output type not recognized. Choose between 'all' and 'preferrred'. Aborting program.");
            exit(1);    
        }
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

double ModuleInterface::getExpectedOccupancy(int queueIndex){
    if (results.evaluated){
        return(results.expectedOccupancy[queueIndex]);
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

double ModuleInterface::getExpOccFraction(int queueIndex){
    if (results.evaluated){
        return(results.expOccFraction[queueIndex]);
    }else{
        py::print("The model has not been evaluated. Aborting program.");
        exit(1);
    }
}

void ModuleInterface::runCalculations(){

    //--------------------------
    //MODEL DATA
    //--------------------------

    //check if input provided
    if(!parametersAvail()){
        py::print("The input parameters are missing. Use importData() to import the parameters. Aborting program.");
        exit(1);
    }
    //indices of queues to be evaluated by the model
    if (settings.evalAll){
        evaluateAllQueues();        
    }

    //check parameters are feasible
    checkParameters();    
    
    //--------------------------
    //MODEL EVALUATION
    //--------------------------

    //control verbose from model (deactivate cout)
    if (!settings.verbose){
        cout.setstate(std::ios_base::failbit);
    }

    //create the model object
    Model mdl(data.arrivalRates,data.serviceTimes,data.capacity,
            data.relocationProbabilities,data.preferredQueue,
            settings.evaluatedQueue,settings.modelType,
            settings.equalizeService);
    
    //additional settings
    if (settings.seed!=-1){
        mdl.setSeed(settings.seed);
    }
    if (settings.burnIn!=-1){
        mdl.setBurnIn(settings.burnIn);
    }
    if (settings.minimumSimulationTime!=-1){
        mdl.setMinimumSimulationTime(settings.minimumSimulationTime);
    }
    if (settings.minSamples!=-1){
        mdl.setMinSamples(settings.minSamples);
    }
    if (settings.hyperStates[0]!=-1 && settings.hyperStates[1]!=-1){
        mdl.setHyperStates(settings.hyperStates[0],settings.hyperStates[1]);
    }
    mdl.setSimTolerance(settings.accTol);
    mdl.setAccuracySampleType(settings.accSampleType);
    
    //now evaluate the model
    mdl.runModel();

    cout.clear(); //reactivate cout
    
    //--------------------------
    // GET THE RESULTS
    //--------------------------

    //queue density distribution    
    results.queueDenDist.resize(data.capacity.size());
    results.queueDenDistPref.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueDenDist[settings.evaluatedQueue[i]].resize(mdl.queueDenDist[settings.evaluatedQueue[i]].size());
        results.queueDenDistPref[settings.evaluatedQueue[i]].resize(mdl.queueDenDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueDenDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueDenDist[settings.evaluatedQueue[i]][j] = mdl.queueDenDist[settings.evaluatedQueue[i]][j];
            results.queueDenDistPref[settings.evaluatedQueue[i]][j] = mdl.queueDenDistPref[settings.evaluatedQueue[i]][j];
        }
    }

    //queue frequency distribution    
    results.queueFreqDist.resize(data.capacity.size());
    results.queueFreqDistPref.resize(data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueFreqDist[settings.evaluatedQueue[i]].resize(mdl.queueFreqDist[settings.evaluatedQueue[i]].size());
        results.queueFreqDistPref[settings.evaluatedQueue[i]].resize(mdl.queueFreqDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueFreqDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueFreqDist[settings.evaluatedQueue[i]][j] = mdl.queueFreqDist[settings.evaluatedQueue[i]][j];
            results.queueFreqDistPref[settings.evaluatedQueue[i]][j] = mdl.queueFreqDistPref[settings.evaluatedQueue[i]][j];
        }
    }
    
    //shortage probabilities
    results.shortageProbability.resize(data.capacity.size(),-1);
    results.shortageProbabilityPref.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.shortageProbability[settings.evaluatedQueue[i]] = mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.shortageProbabilityPref[settings.evaluatedQueue[i]] = mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
    }

    //availability probabilities
    results.availProbability.resize(data.capacity.size(),-1);
    results.availProbabilityPref.resize(data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.availProbability[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.availProbabilityPref[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
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
    
    results.evaluated = true;
}