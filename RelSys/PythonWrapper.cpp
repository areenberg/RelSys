#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>

#include "ModuleInterface.h"

using namespace std;
namespace py = pybind11;


PYBIND11_MODULE(relsys, m) {

    py::class_<ModuleInterface>(m, "model")
        .def(py::init<>())    
        .def("input",&ModuleInterface::importData,"Import the input data to the model.",py::arg("arr"),py::arg("ser"),py::arg("cap"),py::arg("relProb"),py::arg("prefQ")) //import data
        .def("setType",&ModuleInterface::setType,"Set the method to use in the evaluation of the model (auto, simulation, approximation).",py::arg("mdltype")) //model settings
        .def("queuesEval",&ModuleInterface::queuesToEvaluate,"Set the indices of queues to evaluate.",py::arg("qEvalIdx"))
        .def("equalize",&ModuleInterface::equalizeService,"Specify if service times should be equalized and loads correspondingly adjusted (True=On, False=Off).",py::arg("equalize"))
        .def("setVerbose",&ModuleInterface::setVerbose,"Control verbose (True=On, False=Off)",py::arg("set"))
        .def("setSeed",&ModuleInterface::setSeed,"Set the seed.",py::arg("sd"))
        .def("setAccSamType",&ModuleInterface::setAccuracySampleType,"Set the accuracy estimation type for the simulation (preferred, all).",py::arg("stype"))
        .def("setSimTolerance",&ModuleInterface::setSimulationTolerance,"Set the tolerance level for the accuracy estimation in the simulation.",py::arg("tol"))
        .def("setBurnIn",&ModuleInterface::setBurnIn,"Set the burn-in time of the simulation.",py::arg("bin"))
        .def("setSimTime",&ModuleInterface::setMinimumSimulationTime,"Set the simulation time.",py::arg("mnTime"))
        .def("setSamples",&ModuleInterface::setMinSamples,"Set the minimum number of open/shortage samples.",py::arg("mnSamples"))
        .def("setHyperPhases",&ModuleInterface::setHyperStates,"Set the number of phases in the hyper-exponential distributions accounting for the open/shortage time.",py::arg("openStates"),py::arg("blockedStates"))
        .def("run", &ModuleInterface::runCalculations,"Evaluate the model using the input parameters.") //run calculations
        .def("getDensity",&ModuleInterface::getDensityDistribution,"Return the density distribution of a queue.",py::arg("queueIndex")=0,py::arg("type")="preferred") //return results
        .def("getFreq",&ModuleInterface::getFrequencyDistribution,"Return the frequency distribution of a queue.",py::arg("queueIndex")=0,py::arg("type")="preferred")
        .def("getShortageProb",&ModuleInterface::getShortageProbability,"Return the shortage probability of a queue.",py::arg("queueIndex")=0,py::arg("type")="preferred")
        .def("getAvailProb",&ModuleInterface::getAvailProbability,"Return the probability that at least one server is available.",py::arg("queueIndex")=0,py::arg("type")="preferred")
        .def("getExpOccupany",&ModuleInterface::getExpectedOccupancy,"Return the expected number of occupied servers.",py::arg("queueIndex")=0)
        .def("getExpOccFraction",&ModuleInterface::getExpOccFraction,"Return the expected fraction of occupied servers.",py::arg("queueIndex")=0)
        .def("getArrivalRates",&ModuleInterface::getArrivalRates,"Return the imported arrival rates.") //return imported variables
        .def("getServiceTimes",&ModuleInterface::getServiceTimes,"Return the imported service times.")
        .def("getCapacity",&ModuleInterface::getCapacity,"Return the imported capacities.")
        .def("getReloc",&ModuleInterface::getRelocationProbabilities,"Return the imported relocation probabilities.")
        .def("getPreferredQueue",&ModuleInterface::getPreferredQueue,"Return the imported preferred queues.");
        
}