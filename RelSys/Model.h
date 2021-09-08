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

    bool ReadFromFile(const char* fileName);
};