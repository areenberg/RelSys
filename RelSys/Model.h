#pragma once

#include<vector>

struct PatientType
{
    double arrivalRate;
    double meanLengthOfStay;
    std::vector<double> relocationProbabilities;
    int preferredWard;
};

struct Ward
{
    int capacity;
};

struct Model
{
    std::vector<PatientType> patientTypes;
    std::vector<Ward> wards;
    const char* simulationMode;

    int seed;
    int burnInTime;
    int minTime;
    int minSamples;
    int wardIndex;
    int openHyperStates;
    int blockedHyperStates;

    bool ReadFromFile(const char* fileName);
};