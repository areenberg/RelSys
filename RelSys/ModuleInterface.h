#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>

//MODEL TYPES
#include "Model.h"

using namespace std;
namespace py = pybind11;

using namespace std;

#ifndef MODULEINTERFACE_H
#define MODULEINTERFACE_H

class ModuleInterface {
public:

    ModuleInterface();
    ModuleInterface(const ModuleInterface& orig);
    virtual ~ModuleInterface();

    //input data structure
    struct Data{
        vector<double> arrivalRates;
        vector<double> serviceTimes;
        vector<int> capacity;
        vector<vector<double>> relocationProbabilities;
        vector<int> preferredQueue;
    } data;

    //model settings structure
    struct Settings{
        bool verbose=false;
        bool evalAll=true;
        bool equalizeService=false;
        vector<int> evaluatedQueue;
        string modelType="simulation";

        int seed=-1;
        double burnIn=-1;
        double minimumSimulationTime=-1;
        int minSamples=-1;
        vector<int> hyperStates = {-1,-1};
        double accTol=5e-3;
        string accSampleType="preferred";
    } settings;

    //results
    struct Results{
        bool evaluated=false;
        vector<double> shortageProbability,shortageProbabilityPref;
        vector<double> availProbability,availProbabilityPref;
        vector<vector<double>> queueDenDist,queueDenDistPref;
        vector<vector<int>> queueFreqDist,queueFreqDistPref;
        vector<double> expectedOccupancy;
        vector<double> expOccFraction;
    } results;

    void importData(py::list arr, py::list ser, py::list cap,
    py::list relProb, py::list prefQ);
    void setAccuracySampleType(string stype);
    void setSimulationTolerance(double tol);
    void setSeed(int sd);
    void setBurnIn(double bn);
    void setMinimumSimulationTime(double mnTime);
    void setMinSamples(int mnSamples);
    void setHyperStates(int openStates, int blockedStates);
    bool parametersAvail();
    void setVerbose(bool set);
    void checkParameters();
    void evaluateAllQueues();
    void queuesToEvaluate(py::list qEvalIdx);
    void setType(string mdltype);
    void equalizeService(bool equalize);
    py::list getArrivalRates();
    py::list getServiceTimes();
    py::list getCapacity();
    py::list getRelocationProbabilities();
    py::list getPreferredQueue();
    py::list getDensityDistribution(int queueIndex, string type);
    py::list getFrequencyDistribution(int queueIndex, string type);
    double getShortageProbability(int queueIndex, string type);
    double getAvailProbability(int queueIndex, string type);
    double getExpectedOccupancy(int queueIndex);
    double getExpOccFraction(int queueIndex);
    void runCalculations();

private:

};

#endif /* MODULEINTERFACE_H */