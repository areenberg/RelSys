#include "tinyxml.h"
#include "Model.h"


bool Model::ReadFromFile(const char* fileName)
{
    TiXmlDocument doc(fileName);
    if (!doc.LoadFile()) return false;

    TiXmlHandle docHandle(&doc);

    TiXmlElement* l_pPatient = docHandle.FirstChild("Model").FirstChild("PatientTypes").Child("PatientType", 0).ToElement();

    while (l_pPatient)
    {
        PatientType patient;

        // Arrival rate
        {
            TiXmlElement* l_pArrivalRate = l_pPatient->FirstChildElement("ArrivalRate");
            if (NULL != l_pArrivalRate)
            {
                patient.arrivalRate = strtod(l_pArrivalRate->GetText(), NULL);
            }
        }

        // Length of stay
        {
            TiXmlElement* l_pLengthOfStay = l_pPatient->FirstChildElement("LengthOfStay");
            if (NULL != l_pLengthOfStay)
            {
                patient.meanLengthOfStay = strtod(l_pLengthOfStay->GetText(), NULL);
            }
        }

        // Preferred ward
        {
            TiXmlElement* l_pPreferredWard = l_pPatient->FirstChildElement("PreferredWard");
            if (NULL != l_pPreferredWard)
            {
                patient.preferredWard = strtod(l_pPreferredWard->GetText(), NULL);
            }
        }

        // Relocation probabilities
        {
            TiXmlElement* l_pRelocationProbability = l_pPatient->FirstChild("RelocationProbabilities")->FirstChild("Probability")->ToElement();
            while (l_pRelocationProbability)
            {
                if (NULL != l_pRelocationProbability)
                {
                    double probability = strtod(l_pRelocationProbability->GetText(), NULL);
                    patient.relocationProbabilities.push_back(probability);
                }

                l_pRelocationProbability = l_pRelocationProbability->NextSiblingElement();
            }

        }

        patientTypes.push_back(patient);

        l_pPatient = l_pPatient->NextSiblingElement("PatientType");
    }

    TiXmlElement* l_pWard = docHandle.FirstChild("Model").FirstChild("Wards").Child("Ward", 0).ToElement();

    while (l_pWard)
    {
        Ward ward;

        // Capacity
        {
            TiXmlElement* l_pCapacity = l_pWard->FirstChildElement("Capacity");
            if (NULL != l_pCapacity)
            {
                ward.capacity = atoi(l_pCapacity->GetText());
            }
        }

        wards.push_back(ward);

        l_pWard = l_pWard->NextSiblingElement("Ward");
    }

    TiXmlElement* l_pSimulationMode = docHandle.FirstChild("Model").FirstChild("SimulationMode").ToElement();
    simulationMode = l_pSimulationMode->GetText();

    return true;
}