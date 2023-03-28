
/*
 * File:   main.cpp
 * Author: Anders Reenberg Andersen
 *
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Model.h"

using namespace std;


//--------------------------
//STRUCTURES
//--------------------------

//input data structure (including example values)
struct Data {
    vector<double> arrivalRates = {0.8,2.5,0.6,2.8};
    vector<double> serviceTimes {10,5,10,8};
    vector<int> capacity = {15,20,10,30};
    vector<vector<double>> relocationProbabilities = {{0.0,0.4,0.1,0.5},
                                                      {0.3,0.0,0.5,0.0},
                                                      {0.0,0.5,0.0,0.5},
                                                      {0.2,0.3,0.5,0.0}};
    vector<int> preferredQueue = {0,1,2,3};
} data;

//model settings structure
struct Settings {
    bool verbose=false;
    bool evalAll=true;
    bool equalizeService=false;
    vector<int> evaluatedQueue;
    string modelType="simulation";

    int seed=-1;
    double burnIn=-1;
    double minimumSimulationTime=-1;
    int minSamples=-1;
    vector<int> hyperStates = {2,1};
    double accTol=5e-3;
    string accSampleType="preferred";
} settings;

//results
struct Results {
    bool evaluated=false;
    vector<double> shortageProbability,shortageProbabilityPref;
    vector<double> availProbability,availProbabilityPref;
    vector<vector<double>> queueDenDist,queueDenDistPref;
    vector<vector<int>> queueFreqDist,queueFreqDistPref;
    vector<double> expectedOccupancy;
    vector<double> expOccFraction;
} results;

//--------------------------
//METHODS
//--------------------------


//ACTION METHODS

void setVerbose(bool set){
    settings.verbose=set;
}

void evaluateAllQueues(){
    settings.evaluatedQueue.resize(::data.capacity.size());
    for (int i=0; i<::data.capacity.size(); i++){
        settings.evaluatedQueue[i] = i;
    }
}

void queuesToEvaluate(vector<int> qEvalIdx){
    //set the indices of queues to
    //evaluate
    settings.evalAll=false;
    settings.evaluatedQueue = qEvalIdx;
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

void setAccuracySampleType(string stype){
    settings.accSampleType = stype;    
}

void setSimulationTolerance(double tol){
    settings.accTol = tol;
}

void setSeed(int sd){
    settings.seed = sd;
}

void setBurnIn(double bn){
    if (bn>0){
        settings.burnIn=bn;
    }else{
        cout << "Burn-in time must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setMinimumSimulationTime(double mnTime){
    if (mnTime>0){
        settings.minimumSimulationTime=mnTime;
    }else{
        cout << "Simulation time must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setMinSamples(int mnSamples){
    if (mnSamples>0){
        settings.minSamples=mnSamples;
    }else{
        cout << "Min. number of open/shortage samples must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setOpenHyperStates(int openStates){
    if (openStates>0){
        settings.hyperStates[0] = openStates;
    }else{
        cout << "The number of open and blocked states must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void setBlockedHyperStates(int blockedStates){
    if (blockedStates>0){
        settings.hyperStates[1] = blockedStates;
    }else{
        cout << "The number of open and blocked states must be larger than 0. Aborting program." << endl;
        exit(1);
    }
}

void printHelp(){
    cout << "Usage: relsys.exe [options]" << endl <<
    endl << "A demo is loaded if no options are provided." << endl << endl;
    cout << "Options include:" << endl << "-v   Activate verbose." << endl;
    cout << "-arr <filename>   A space-separated vector containing the arrival rates of each customer type." << endl;
    cout << "-ser <filename>   A space-separated vector containing the service times of each customer type." << endl;
    cout << "-cap <filename>   A space-separated vector containing the capacities of each queue." << endl;
    cout << "-rel <filename>   A space-separated matrix. Each element denotes the probability that customer type i will try to use queue j when the preferred queue is in shortage." << endl;
    cout << "-prq <filename>   A space-separated vector containing the indices of the queues that are preferred by each customer type." << endl;
    cout << "-evq <filename>   A space-separated vector containing the indices of the queues that are evaluated by the program." << endl;
    
    cout << "-m <value>   Select the model type (simulation (default), approximation, auto)." << endl;
    cout << "-sd <value>    Set a seed for the model." << endl;
    cout << "-tol <value>    Set the tolerance for the automatic termination of the simulation." << endl;
    cout << "-asamt <value>    Set the accuracy evaluation type for the automatic termination of the simulation (preferred (default), all)." << endl;
    cout << "-bin <value>    Set the burn-in time for the simulation." << endl;
    cout << "-stime <value>    Set the simulation time for the simulation." << endl;
    cout << "-msam <value>    Set the minimum number of open/shortage samples for the approximation." << endl;
    cout << "-eq    Equalize the service-rates." << endl;
    cout << "-ops <value>    The number of open states for the approximation." << endl;
    cout << "-bla <value>    The number of shortage/blocked states for the approximation." << endl;
}

void checkParameters(){
    //checks if parameters are feasible

    if (::data.arrivalRates.size() != ::data.serviceTimes.size()){
        cout << "The number of arrival rates must equal the number of service times. Aborting program." << endl;
        exit(1);
    }
    if (::data.arrivalRates.size() != ::data.relocationProbabilities.size()){
        cout << "The number of arrival rates, and service times, must equal the number of rows in the relocation matrix. Aborting program." << endl;
        exit(1);
    }
    if (::data.arrivalRates.size() != ::data.preferredQueue.size()){
        cout << "The number of arrival rates, and service times, must equal the length of the vector specifying the preferred queues. Aborting program." << endl;
        exit(1);
    }
    if (settings.evaluatedQueue.size()>::data.capacity.size() || settings.evaluatedQueue.empty()){
        cout << "The length of the vector of evaluated queues must be between 1 and n_queues. Aborting program." << endl;
        exit(1);
    }        
    for (int i=0; i<::data.relocationProbabilities.size(); i++){
        if (::data.capacity.size() != ::data.relocationProbabilities[i].size()){
            cout << "The length of the capacity vector must equal the number of columns in the relocation matrix. Aborting program." << endl;
            exit(1);
        }
    }        
    for (int i=0; i<::data.arrivalRates.size(); i++){
        if (::data.arrivalRates[i]<=0.0 || ::data.serviceTimes[i]<=0.0 ){
            cout << "Arrival rates and service times must be larger than 0.0. Aborting program." << endl;
            exit(1);
        }
        for (int j=0; j<::data.relocationProbabilities[i][j]; j++){
            if (::data.relocationProbabilities[i][j]<0.0 || ::data.relocationProbabilities[i][j]>1.0){
                cout << "Values in the relocation matrix must be between 0.0 and 1.0. Aborting program." << endl;
                exit(1);
            }
        }
    }
    for (int i=0; i<::data.preferredQueue.size(); i++){
        if (::data.preferredQueue[i]<0 || ::data.preferredQueue[i]>(::data.capacity.size()-1)){
            cout << "Indices in the vector of preferred queues must be between 0 and n_queues-1. Aborting program." << endl;
            exit(1);
        }
    }
    for (int i=0; i<::data.capacity.size(); i++){
        if (::data.capacity[i]<1){
            cout << "The capacity of each queue must be equal to or larger than 1. Aborting program." << endl;
            exit(1);
        }
    }          
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        if (settings.evaluatedQueue[i]<0 || settings.evaluatedQueue[i]>(::data.capacity.size()-1)){
            cout << "Indices in the vector of evaluated queues must lie in the interval between 0 and n_queues-1. Aborting program." << endl;
            exit(1);
        }    
    }
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        for (int j=0; j<settings.evaluatedQueue.size(); j++){
            if (i!=j && settings.evaluatedQueue[i]==settings.evaluatedQueue[j]){
                cout << "All indices in the vector of evaluated queues must be unique. The same queue cannot be evaluated twice. Aborting program." << endl;
                exit(1);
            }
        }    
    }
    double sm;
    for (int i=0; i<::data.relocationProbabilities.size(); i++){
        sm=0;
        for (int j=0; j<::data.relocationProbabilities[i].size(); j++){
            sm+=::data.relocationProbabilities[i][j];
        }
        if (sm>1.0){
            string out = "The sum of the relocation probabilities in row "+to_string(i)+" is equal to "+to_string(sm)+". The sum must be equal to or smaller than 1.0. Aborting program.";
            cout << out << endl;
            exit(1);
        }
    }
    if ((settings.minimumSimulationTime!=-1 && settings.burnIn==-1) || (settings.minimumSimulationTime==-1 && settings.burnIn!=-1)){
        cout << "It is not possible to set the overall simulation time without setting the burn-in time and vice versa. Aborting program." << endl;
        exit(1);
    }
    if (settings.minimumSimulationTime!=-1 && settings.burnIn!=-1 && settings.minimumSimulationTime<=settings.burnIn){
        cout << "The simulation time has to be longer than the burn-in time. Aborting program." << endl;
        exit(1);
    }
    if (settings.accTol<=0.0 || settings.accTol>=1.0){
        cout << "The tolerance level for the accuracy estimation procedure must be between 0 and 1 (e.g. 1e-2). Aborting program." << endl;
        exit(1);
    }        

}

//READ/WRITE FILES
vector<vector<double>> readMatrixFromFile(string fileName) {
    vector<vector<double>> mat;
    ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        cout << "The file " << fileName << " was not found or could not be opened. Aborting program." << endl;
        exit(1);
    }

    int nRows, nCols;
    vector<int> rowLengths;
    string line;

    while (getline(inputFile, line)) {
        int count = 0;
        for (char c : line) {
            if (c == ' ') {
                count++;
            }
        }
        rowLengths.push_back(count + 1);
    }

    nRows = rowLengths.size();
    nCols = rowLengths[0];

     for (int i = 1; i < nRows; i++) {
        if (rowLengths[i] != nCols) {
            cout << "The matrix in file " << fileName << " is not rectangular. Aborting program." << endl;
            exit(1);
        }
    }
    inputFile.close();

    mat.resize(nRows, vector<double>(nCols));
    string val;

    ifstream inputFile2(fileName);
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            inputFile2 >> val;
            mat[i][j] = stod(val);
        }
    }

    inputFile2.close();

    return(mat);
}



// READ ARGUMENTS
double argDouble(string flag, int &argc, char** argv){
    int i=0;
    string str;
    double val;
    do{
        str=argv[i];
        if (str.compare(flag)==0){
            if (i<(argc-1)){
                return(atof(argv[i+1]));
            }else{
                cout << "Missing input value for flag " << flag << endl;
            }       
        }
        i++;
    }while(i<argc && str.compare(flag)!=0);
    return(-1);
}

double argInteger(string flag, int &argc, char** argv){
    int i=0;
    string str;
    int val;
    do{
        str=argv[i];
        if (str.compare(flag)==0){
            if (i<(argc-1)){
                return(stoi(argv[i+1]));
            }else{
                cout << "Missing input value for flag " << flag << endl;
            }       
        }
        i++;
    }while(i<argc && str.compare(flag)!=0);
    return(-1);
}

string argString(string flag, int &argc, char** argv){
    int i=0;
    string str1,str2;
    do{
        str1=argv[i];
        if (str1.compare(flag)==0){
            if (i<(argc-1)){
                str2 = argv[i+1];
                return(str2);
            }else{
                cout << "Missing input value for flag " << flag << endl;
            }       
        }
        i++;
    }while(i<argc && str1.compare(flag)!=0);
    return("NA");
}

bool argActivate(string flag, int &argc, char** argv){
    int i=0;
    string str;
    do{
        str=argv[i];
        if (str.compare(flag)==0){
            return(true);       
        }
        i++;
    }while(i<argc && str.compare(flag)!=0);
    return(false);
}

void argArrivalRates(int &argc, char** argv){
    string fileName = argString("-arr",argc,argv);
    if (fileName.compare("NA")!=0){
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        ::data.arrivalRates.resize(mat[0].size());
        for (int j=0; j<mat[0].size(); j++){
            ::data.arrivalRates[j] = mat[0][j];
        }
    }
}

void argServiceTimes(int &argc, char** argv){
    string fileName = argString("-ser",argc,argv);
    if (fileName.compare("NA")!=0){
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        ::data.serviceTimes.resize(mat[0].size());
        for (int j=0; j<mat[0].size(); j++){
            ::data.serviceTimes[j] = mat[0][j];
        }
    }
}

void argCapacity(int &argc, char** argv){
    string fileName = argString("-cap",argc,argv);
    if (fileName.compare("NA")!=0){
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        ::data.capacity.resize(mat[0].size());
        for (int j=0; j<mat[0].size(); j++){
            ::data.capacity[j] = (int)mat[0][j];
        }
    }
}

void argPreferred(int &argc, char** argv){
    string fileName = argString("-prq",argc,argv);
    if (fileName.compare("NA")!=0){
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        ::data.preferredQueue.resize(mat[0].size());
        for (int j=0; j<mat[0].size(); j++){
            ::data.preferredQueue[j] = (int)mat[0][j];
        }
    }
}

void argEvalQueues(int &argc, char** argv){
    string fileName = argString("-evq",argc,argv);
    if (fileName.compare("NA")!=0){
        settings.evalAll=false;
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        settings.evaluatedQueue.resize(mat[0].size());
        for (int j=0; j<mat[0].size(); j++){
            settings.evaluatedQueue[j] = (int)mat[0][j];
        }
    }
}

void argRelocProbs(int &argc, char** argv){
    string fileName = argString("-rel",argc,argv);
    if (fileName.compare("NA")!=0){
        vector<vector<double>> mat = readMatrixFromFile(fileName);
        ::data.relocationProbabilities.resize(mat.size());
        for (int i=0; i<mat.size(); i++){
            ::data.relocationProbabilities[i].resize(mat[i].size());
            for (int j=0; j<mat[i].size(); j++){
                ::data.relocationProbabilities[i][j] = mat[i][j];    
            }
        }
    }
}

void argModel(int &argc, char** argv){
    string str = argString("-m",argc,argv);
    if (str.compare("NA")!=0){
        setType(str);
    }    
}

void argSeed(int &argc, char** argv){
    int val = argInteger("-sd",argc,argv);
    if (val!=-1){
        setSeed(val);
    }
}

void argOpenStates(int &argc, char** argv){
    int val = argInteger("-ops",argc,argv);
    if (val!=-1){
        setOpenHyperStates(val);
    }
}

void argBlockedStates(int &argc, char** argv){
    int val = argInteger("-bls",argc,argv);
    if (val!=-1){
        setBlockedHyperStates(val);
    }
}

void argAccType(int &argc, char** argv){
    string str = argString("-asamt",argc,argv);
    if (str.compare("NA")!=0){
        setAccuracySampleType(str);
    }    
}

void argTol(int &argc, char** argv){
    double val = argDouble("-tol",argc,argv);
    if (val!=-1){
        setSimulationTolerance(val);
    }    
}

void argBurnIn(int &argc, char** argv){
    double val = argDouble("-bin",argc,argv);
    if (val!=-1){
        setBurnIn(val);
    }    
}

void argSimTime(int &argc, char** argv){
    double val = argDouble("-stime",argc,argv);
    if (val!=-1){
        setMinimumSimulationTime(val);
    }    
}

void argMinSamples(int &argc, char** argv){
    int val = argInteger("-msam",argc,argv);
    if (val!=-1){
        setMinSamples(val);
    }    
}

void argVerbose(int &argc, char** argv){
    bool act = argActivate("-v",argc,argv);
    if (act){
        setVerbose(true);
    }
}

void argEqualize(int &argc, char** argv){
    bool act = argActivate("-eq",argc,argv);
    if (act){
        equalizeService(true);
    }
}

void readArguments(int &argc, char** argv){
    argArrivalRates(argc,argv);
    argServiceTimes(argc,argv);
    argRelocProbs(argc,argv);
    argPreferred(argc,argv);
    argEvalQueues(argc,argv);
    argVerbose(argc,argv);        
    argModel(argc,argv);
    argEqualize(argc,argv);
    argAccType(argc,argv);
    argTol(argc,argv);
    argSeed(argc,argv);
    argBurnIn(argc,argv);
    argSimTime(argc,argv);
    argMinSamples(argc,argv);
    argOpenStates(argc,argv);
    argBlockedStates(argc,argv);
}

void runCalculations(){

    //--------------------------
    //MODEL DATA
    //--------------------------

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
    Model mdl(::data.arrivalRates,::data.serviceTimes,::data.capacity,
            ::data.relocationProbabilities,::data.preferredQueue,
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
    if (settings.hyperStates[0]!=2 || settings.hyperStates[1]!=1){
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
    results.queueDenDist.resize(::data.capacity.size());
    results.queueDenDistPref.resize(::data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueDenDist[settings.evaluatedQueue[i]].resize(mdl.queueDenDist[settings.evaluatedQueue[i]].size());
        results.queueDenDistPref[settings.evaluatedQueue[i]].resize(mdl.queueDenDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueDenDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueDenDist[settings.evaluatedQueue[i]][j] = mdl.queueDenDist[settings.evaluatedQueue[i]][j];
            results.queueDenDistPref[settings.evaluatedQueue[i]][j] = mdl.queueDenDistPref[settings.evaluatedQueue[i]][j];
        }
    }

    //queue frequency distribution    
    results.queueFreqDist.resize(::data.capacity.size());
    results.queueFreqDistPref.resize(::data.capacity.size());
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.queueFreqDist[settings.evaluatedQueue[i]].resize(mdl.queueFreqDist[settings.evaluatedQueue[i]].size());
        results.queueFreqDistPref[settings.evaluatedQueue[i]].resize(mdl.queueFreqDistPref[settings.evaluatedQueue[i]].size());
        for (int j=0; j<mdl.queueFreqDist[settings.evaluatedQueue[i]].size(); j++){
            results.queueFreqDist[settings.evaluatedQueue[i]][j] = mdl.queueFreqDist[settings.evaluatedQueue[i]][j];
            results.queueFreqDistPref[settings.evaluatedQueue[i]][j] = mdl.queueFreqDistPref[settings.evaluatedQueue[i]][j];
        }
    }
    
    //shortage probabilities
    results.shortageProbability.resize(::data.capacity.size(),-1);
    results.shortageProbabilityPref.resize(::data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.shortageProbability[settings.evaluatedQueue[i]] = mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.shortageProbabilityPref[settings.evaluatedQueue[i]] = mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
    }

    //availability probabilities
    results.availProbability.resize(::data.capacity.size(),-1);
    results.availProbabilityPref.resize(::data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.availProbability[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbability[settings.evaluatedQueue[i]];
        results.availProbabilityPref[settings.evaluatedQueue[i]] = 1.0-mdl.blockingProbabilityPref[settings.evaluatedQueue[i]];
    }

    //expected occupancy
    results.expectedOccupancy.resize(::data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expectedOccupancy[settings.evaluatedQueue[i]] = mdl.expectedOccupancy[settings.evaluatedQueue[i]];
    }

    //expected fraction of customers occupied
    results.expOccFraction.resize(::data.capacity.size(),-1);
    for (int i=0; i<settings.evaluatedQueue.size(); i++){
        results.expOccFraction[settings.evaluatedQueue[i]] = mdl.expOccFraction[settings.evaluatedQueue[i]];
    }
    
    results.evaluated = true;

}

int main(int argc, char** argv) {

    //--------------------------
    // READ ARGUMENTS
    //--------------------------
    
    readArguments(argc,argv);
    
    if (argActivate("-help",argc,argv)){
        printHelp();
    }else{

        //--------------------------
        // EVALUATE MODEL
        //--------------------------
    
        runCalculations();

        //--------------------------
        // DISPLAY OR SAVE RESULTS
        //--------------------------
    
        //CODE HERE
    
    }

    return 0;
}