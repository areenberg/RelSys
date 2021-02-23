/*
 * Copyright 2020 Anders Reenberg Andersen.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* 
 * File:   RelocSimulation.cpp
 * Author: Anders Reenberg Andersen
 * 
 * Created on November 13, 2020, 2:02 PM
 */

#include "RelocSimulation.h"
#include "Patient.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <limits>
#include <random>

using namespace std;

RelocSimulation::RelocSimulation(int nW, WardData * wards):
wards_pointer(wards),
nWards(nW)        
{
    initializeSystem();
}

RelocSimulation::RelocSimulation(const RelocSimulation& orig) {
}

RelocSimulation::~RelocSimulation() {
}

void RelocSimulation::initializeSystem(){
    calculateArrivalRates();
    
    //setup for the random number generator
    random_device rd;
    simSeed = rd();
    mt19937 rgen(simSeed);
    uniform_real_distribution<> dis(0,1);
}

void RelocSimulation::calculateArrivalRates(){
    //calculates the arrival rate into each bin of the system
    
    arrivalRateMatrix.resize(nWards);
    for (int widx=0; widx<nWards; widx++){
        arrivalRateMatrix[widx].resize(nWards,0);
        for (int pidx=0; pidx<nWards; pidx++){
            if (widx==pidx){
                arrivalRateMatrix[widx][pidx] = getWardArrivalRate(widx);
            }else{
                arrivalRateMatrix[widx][pidx] = getWardArrivalRate(pidx)*
                        getWardRelocationProbabilities(pidx)[widx];
            }
        }
    }
    
}

void RelocSimulation::setSeed(int seed){
    //srand(seed);
    simSeed = seed;
    mt19937 rgen(simSeed);
}

void RelocSimulation::initializeArrivalTimes(double currentClock){
    //allocate memory and fill in the entire nextArrivalTime matrix
    
    nextArrivalTime.clear();
    nextArrivalTime.resize(nWards);
    for (int i=0; i<nWards; i++){
        nextArrivalTime[i].resize(nWards,currentClock);
    }
    
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=0; pidx<nWards; pidx++){
            nextArrivalTime[widx][pidx] += randomExponential(arrivalRateMatrix[widx][pidx]);
        }
    }
    
}

void RelocSimulation::updateArrivalTime(int widx, int pidx){
    //derives the next arrival time for a single ward-patient pair
    
    nextArrivalTime[widx][pidx] += randomExponential(arrivalRateMatrix[widx][pidx]);
    
}

void RelocSimulation::generateArrivalList(int length, double currentClock){
    //generate a long list of future arrivals with arrival times offset
    //by currentClock
    
    initializeArrivalTimes(currentClock);
    
    //append patients to the list
    double arrClock, serTime, mn;
    int min_widx, min_pidx;
    //vector<int> genWidx(nWards,0);
    //vector<int> genPidx(nWards,0);
    
    for (int i=0; i<length; i++){
        
        mn = numeric_limits<double>::max();
        
        //find minimum arrival time        
        for (int widx=0; widx<nWards; widx++){
            for (int pidx=0; pidx<nWards; pidx++){
                if (nextArrivalTime[widx][pidx]<mn){
                    mn = nextArrivalTime[widx][pidx];
                    min_widx = widx; min_pidx = pidx;
                }
            }
        }        
        
        
        //append patient
        arrClock = nextArrivalTime[min_widx][min_pidx];
        serTime = randomExponential(getWardServiceRate(min_pidx)); //<<-- generates the random service time
        //serTime = randomLogNormal((1/getWardServiceRate(min_pidx)),(1/getWardServiceRate(min_pidx))); //<<-- change to log-normal service time
        arrival_array[i] = Patient(arrClock,serTime,min_widx,min_pidx);
        
        //update ward-patient
        updateArrivalTime(min_widx,min_pidx);
    }
    
    
}

void RelocSimulation::simulate(double bIn, double minTime, int minSamples){
    //burnIn is the burn-in time.
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
    //variables
    double nextArrClock, nextSerClock;
    int patientArraySize, maxOcc, inService, arrIdx, 
            serIdx, targetWard;
    bool succeeded, note = false;
    //bool adm;
    
    vector<vector<int>> wardOccupancy;
    vector<int> nOpenTimeSamples(nWards,0);
    vector<int> nBlockedTimeSamples(nWards,0);
    wardStateClocks.resize(nWards,0);
    
    //--------------------
    //Initialize
    //--------------------
    
    //burn-in time
    burnIn = bIn;
    
    //simulation clock
    clock = 0;
    
    //patient arrivals
    patientArraySize = 1e5;
    arrival_array = new Patient[patientArraySize];
    generateArrivalList(patientArraySize,clock);
    arrIdx = 0;
    
    //patients in service
    maxOcc = 0; inService = 0;
    for (int widx=0; widx<nWards; widx++){
        maxOcc += getWardCapacity(widx);
    }
    service_array = new Patient[maxOcc];
    
    //occupancy of the system
    vector<int> capUse(nWards,0);
    
    //tracking
    openTimes.resize(nWards);
    blockedTimes.resize(nWards);
    //initialize the frequency/density
    //ward-patient distributions.
    initFreqDenDist();  
    
    //-------------------
    //Simulate
    //-------------------
    cout << "Running simulation..." << flush;
    while (arrIdx<patientArraySize && (clock<minTime || 
            minTimeSamples()<minSamples)){
        
        //clock at next arrival
        nextArrClock = arrival_array[arrIdx].arrivalClock;
        
//        cout << "sim. clock = " << clock << endl;
//        for (int i=0; i<nWards; i++){
//            cout << wardStateClocks[i] << ", ";
//        }
//        cout << endl;
//        cout << "admission widx=" << arrival_array[arrIdx].wardTarget << ", pidx=" << arrival_array[arrIdx].patientType << ", adm.clock=" << nextArrClock << endl; 
//        cout << "cap use=" << capUse[arrival_array[arrIdx].wardTarget] << ", cap lim=" << getWardCapacity(arrival_array[arrIdx].wardTarget) << endl;
//    
//        cout << "old array:" << endl;
//        for (int i=0; i<inService; i++){
//            cout << service_array[i].serviceClock << "," << service_array[i].wardTarget << "," << service_array[i].patientType << endl;    
//        }
    
        
        if (inService>0){
            //clock at next discharge
            serIdx = nextServiceIdx(inService);
            nextSerClock = service_array[serIdx].serviceClock;
//            cout << "inService=" << inService << ", serIdx=" << serIdx << ", nextSerClock=" << nextSerClock << endl;
            
            if (nextSerClock<nextArrClock){
                clock = nextSerClock;
                targetWard = service_array[serIdx].wardTarget;
                succeeded = attemptDischarge(serIdx,inService,capUse);
                if (succeeded){
                    blockedTimeTracking(nBlockedTimeSamples,targetWard,capUse);
                }
            }else{
                clock = nextArrClock;
                succeeded = attemptAdmission(arrIdx,capUse,wardOccupancy,inService);
                if (succeeded){
                    openTimeTracking(nOpenTimeSamples,
                        service_array[inService-1].wardTarget,capUse);
                }
            }
        }else{
            clock = nextArrClock; 
            succeeded = attemptAdmission(arrIdx,capUse,wardOccupancy,inService);
            if (succeeded){
                openTimeTracking(nOpenTimeSamples,
                    service_array[inService-1].wardTarget,capUse);
            }
        }
        
        updateOccupancy(capUse,wardOccupancy,inService);
        
//        cout << "new array:" << endl;
//        for (int i=0; i<inService; i++){
//            cout << service_array[i].serviceClock << "," << service_array[i].wardTarget << "," << service_array[i].patientType << endl; 
//        }
//        cout << "----------" << endl;
    
        
        //generate a new batch of patient arrivals
        if (arrIdx==(patientArraySize-1)){
            generateArrivalList(patientArraySize,clock);
            arrIdx = 0;
        }else if(note==false && clock>minTime && minTimeSamples()<minSamples){
            cout << "still collecting samples..." << flush;
            note = true;
        }else if(note==false && clock<minTime && minTimeSamples()>=minSamples){
            cout << "state-time samples collected..." << flush;
            note = true;
        }
        
    }
    
    freqToDensity();
    
    cout << "done." << endl;
    
    //re-sample time tracking
    subsetTimeSamples(minSamples);
    
    //print open/blocked time sample sizes
    printTimeSamples();
    
}

void RelocSimulation::subsetTimeSamples(int &minSamples){
    
    vector<int> idx; //indices to remove or move
    vector<double> x(minSamples,0);
    int from=0,to;
    
    //open time samples
    for (int i=0; i<openTimes.size(); i++){
        if(openTimes[i].size()>minSamples){
            x.clear(); x.resize(minSamples,0);
            idx.clear();
            to=openTimes[i].size()-1;
            idx = randomIndices(from,to,minSamples);
            for (int j=0; j<x.size(); j++){
                x[j] = openTimes[i][idx[j]];
            }
            openTimes[i].resize(minSamples,0);
            for (int j=0; j<minSamples; j++){
                openTimes[i][j] = x[j];
            }
                
        }    
    }
    
    //blocked time samples
    for (int i=0; i<blockedTimes.size(); i++){
        if(blockedTimes[i].size()>minSamples){
            x.clear(); x.resize(minSamples,0);
            idx.clear();
            to = blockedTimes[i].size()-1;
            idx = randomIndices(from,to,minSamples);
            for (int j=0; j<x.size(); j++){
                x[j] = blockedTimes[i][idx[j]];
            }
            blockedTimes[i].resize(minSamples,0);
            for (int j=0; j<minSamples; j++){
                blockedTimes[i][j] = x[j];
            }
                
        }    
    }
    
}

vector<int> RelocSimulation::randomIndices(int &from, int &to, int &len){
    
    int l,added,nn,
            span,idx;
    bool rerun;
    
    span = (to-from)+1;
    l = min(len,span);
    vector<int> vec(l,0);
    if (l<len || len==span){
        for (int i=0; i<vec.size(); i++){
            vec[i] = i;
        }
    }else{
        added=0;
        while (added<l){
            rerun = true;
            while (rerun){
                nn = floor(randomUniform()*span+from);
                if (added==0){
                    vec[added] = nn;
                    added++;
                    rerun=false;
                }else{
                    idx=0;
                    while (idx<l && vec[idx]!=nn){
                        idx++;
                    }
                    if (idx==l){
                        vec[added] = nn;
                        added++;
                        rerun=false;
                    }
                }
            }
        }
    }
    
    return(vec);
}


void RelocSimulation::printTimeSamples(){
    
    cout << "----------------------------\nSIM STATS:" << endl;
    cout << "*Number of open-time samples*" << endl;
    for (int i=0; i<openTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << openTimes[i].size();
        if (i<(openTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << "\n" << endl;
    cout << "*Number of blocked-time samples*" << endl;
    for (int i=0; i<blockedTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << blockedTimes[i].size();
        if (i<(blockedTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << endl;
    cout << "----------------------------\n" << endl;
    
}

void RelocSimulation::initFreqDenDist(){
    
    wardFreqDist.resize(nWards);//marginal frequency distribution
    freqDist.resize(nWards); //frequency distribution
    denDist.resize(nWards); //density distribution
    
    for (int widx=0; widx<nWards; widx++){
        int c = getWardCapacity(widx)+1;
        wardFreqDist[widx].resize(c,0);
        freqDist[widx].resize(nWards);
        denDist[widx].resize(nWards);
        for (int pidx=0; pidx<nWards; pidx++){
            freqDist[widx][pidx].resize(c,0);
            denDist[widx][pidx].resize(c,0);
        }
    }
    
}

void RelocSimulation::updateOccupancy(vector<int> &capUse, vector<vector<int>> &wardOccupancy, int &inService){
    
    //ward occupancy (all patients)
    capUse.clear();
    capUse.resize(nWards,0);
    
    for (int i=0; i<inService; i++){
        capUse[service_array[i].wardTarget]++;
    }
    
    //ward occupancy (each patient type)
    wardOccupancy.clear();
    wardOccupancy.resize(nWards);
    for (int widx=0; widx<nWards; widx++){
        wardOccupancy[widx].resize(nWards,0);
    }
    for (int i=0; i<inService; i++){
        wardOccupancy[service_array[i].wardTarget][service_array[i].patientType]++;
    }
    
}

void RelocSimulation::occupancyDistTracking(vector<vector<int>> &wardOccupancy,
    vector<int> &capUse, int &targetWard, int &patientType){
    
    if (clock>burnIn){
        if ( (targetWard!=patientType && getWardCapacity(patientType)==capUse[patientType]) ||
                targetWard==patientType ){
            //record distribution over all patient types in the ward
            int wardSum=0;
            for (int pidx=0; pidx<nWards; pidx++){
                freqDist[targetWard][pidx][wardOccupancy[targetWard][pidx]]++;
                wardSum += wardOccupancy[targetWard][pidx];
            }
            wardFreqDist[targetWard][wardSum]++;
        }
    }
}

void RelocSimulation::freqToDensity(){
    
    vector<vector<int>> sm(nWards);
    
    for (int widx=0; widx<nWards; widx++){
        sm[widx].resize(nWards,0);
        for (int pidx=0; pidx<nWards; pidx++){
            for (int i=0; i<freqDist[widx][pidx].size(); i++){
                sm[widx][pidx] += freqDist[widx][pidx][i];
            }
        }
    }
    
    for (int widx=0; widx<nWards; widx++){
        for (int pidx=0; pidx<nWards; pidx++){
            for (int i=0; i<denDist[widx][pidx].size(); i++){
                denDist[widx][pidx][i] = (double) freqDist[widx][pidx][i]/ (double) sm[widx][pidx];
            }
        }
    }
    
}

void RelocSimulation::openTimeTracking(
        vector<int> &nOpenTimeSamples, int &targetWard, vector<int> &capUse){
    
    if ((capUse[targetWard]+1)==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
            openTimes[targetWard].push_back(sample);
            nOpenTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}
    
void RelocSimulation::blockedTimeTracking(
        vector<int> &nBlockedTimeSamples, int &targetWard, vector<int> &capUse){
    
    if (capUse[targetWard]==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
//            if (sample==0){
//                cout << "SAMPLE 0" << endl;
//            }
            blockedTimes[targetWard].push_back(sample);
            nBlockedTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}



bool RelocSimulation::attemptDischarge(int &serIdx, int &inService, vector<int> &capUse){
    
    //subtract patient from service
    inService--;
    
    updateServiceArray(serIdx,inService);
    
    //cout << "attempted discharge success" << endl;
    
    return(true);
}

bool RelocSimulation::attemptAdmission(int &arrIdx, vector<int> &capUse, 
        vector<vector<int>> &wardOccupancy, int &inService){
    
    bool Ok = false;
    if (arrival_array[arrIdx].wardTarget==arrival_array[arrIdx].patientType &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }else if(arrival_array[arrIdx].wardTarget!=arrival_array[arrIdx].patientType &&
            capUse[arrival_array[arrIdx].patientType]==getWardCapacity(arrival_array[arrIdx].patientType) &&
            capUse[arrival_array[arrIdx].wardTarget]<getWardCapacity(arrival_array[arrIdx].wardTarget)){
        Ok = true;
    }
    
    if (Ok){
        
        //insert into service array
        service_array[inService] = Patient(arrival_array[arrIdx].arrivalClock,arrival_array[arrIdx].serviceTime,
                arrival_array[arrIdx].wardTarget,arrival_array[arrIdx].patientType);
        //calculate and insert clock at discharge
        service_array[inService].serviceClock = clock + service_array[inService].serviceTime;
        
        //adjust number of patients currently in service
        inService++;
    }
    
    //track regardless of admission accepted or rejected
    occupancyDistTracking(wardOccupancy,capUse,arrival_array[arrIdx].wardTarget,
                    arrival_array[arrIdx].patientType);
    
    //move to next arrival
    arrIdx++;
    
    return(Ok);
}
    
void RelocSimulation::updateServiceArray(int idx, int &inService){
    //removes a patient in service and adjusts the entire list
    //assumes the patient <<has been subtracted>> from inService
    
    if (idx<inService){
        double arrClock, serTime, serClock;
        int widx, pidx;
        
        for (int i=idx; i<inService; i++){
            arrClock = service_array[i+1].arrivalClock;
            serTime = service_array[i+1].serviceTime;
            serClock = service_array[i+1].serviceClock;
            widx = service_array[i+1].wardTarget;
            pidx = service_array[i+1].patientType;
            
            service_array[i] = Patient(arrClock,serTime,widx,pidx);
            service_array[i].serviceClock = serClock;
        }
    }
    
}

int RelocSimulation::nextServiceIdx(int &inService){
    
    if (inService>1){
        int idx;
        double mn = numeric_limits<double>::max();
        for (int i=0; i<inService; i++){
            if (service_array[i].serviceClock<mn){
               mn = service_array[i].serviceClock;
               idx = i;
            }
        }
//        cout << "serIdx=" << idx << ", clock=" << mn << endl;
        return(idx);
    }else{
        return(0);
    }
    
}

int RelocSimulation::minTimeSamples(){
    
    int mn = numeric_limits<int>::max();
    
    for (int i=0; i<openTimes.size(); i++){
        if (mn>openTimes[i].size()){
            mn = openTimes[i].size();
        }
    }
    for (int i=0; i<blockedTimes.size(); i++){
        if (mn>blockedTimes[i].size()){
            mn = blockedTimes[i].size();
        }
    }
    
    return(mn);
}


double RelocSimulation::randomExponential(double rate){
    //generate a random exponential double
    return(log(1-randomUniform())/(-rate));
}

double RelocSimulation::randomLogNormal(double mean, double std){
    //generate a random log-normal double
    double mu = log(pow(mean,2.0)/(sqrt(pow(mean,2.0)+pow(std,2.0))));
    double sigma = sqrt(log(1+(pow(std,2.0)/pow(mean,2.0))));
    return(logNormalCdfInv(randomUniform(),mu,sigma));
}

double RelocSimulation::logNormalCdfInv(double cdf, double mu, double sigma){
  double logx,x;
  
  if ( cdf < 0.0 || 1.0 < cdf ){
    cerr << "\n";
    cerr << "LOG_NORMAL_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit(1);
  }

  logx = normalCdfInv(cdf,mu,sigma);
  x = exp(logx);

  return x;
}

double RelocSimulation::normalCdfInv(double cdf, double mu, double sigma){
    double x,x2;
  
    if (cdf < 0.0 || 1.0 < cdf){
        cerr << "\n";
        cerr << "NORMAL_CDF_INV - Fatal error!\n";
        cerr << "  CDF < 0 or 1 < CDF.\n";
        exit ( 1 );
    }
    x2 = normal01CdfInv(cdf);
    x = mu + sigma * x2;
    return x;
}

double RelocSimulation::normal01CdfInv(double p){
  double a[8] = {
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770e-1,
    2.27238449892691845833e-2,  7.74545014278341407640e-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550e-1,
    1.48103976427480074590e-1,  1.51986665636164571966e-2,
    5.47593808499534494600e-4,  1.05075007164441684324e-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230e-1,
    2.65321895265761230930e-2,  1.24266094738807843860e-3,
    2.71155556874348757815e-5,  2.01033439929228813265e-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690e-1,
    1.36929880922735805310e-1,  1.48753612908506148525e-2,
    7.86869131145613259100e-4,  1.84631831751005468180e-5,
    1.42151175831644588870e-7,  2.04426310338993978564e-15 };
  double q;
  double r;
  const double r8_huge = 1.0E+30;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;

  if ( p <= 0.0 )
  {
    value = - r8_huge;
    return value;
  }

  if ( 1.0 <= p )
  {
    value = r8_huge;
    return value;
  }

  q = p - 0.5;

  if (fabs ( q ) <= split1){
    r = const1 - q * q;
    value = q * r8polyValueHorner ( 7, a, r ) / r8polyValueHorner ( 7, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = r8_huge;
    }
    else
    {
      r = sqrt ( - log ( r ) );

      if ( r <= split2 ){
        r = r - const2;
        value = r8polyValueHorner ( 7, c, r ) / r8polyValueHorner ( 7, d, r );
       }
       else
       {
         r = r - split2;
         value = r8polyValueHorner ( 7, e, r ) / r8polyValueHorner ( 7, f, r );
      }
    }

    if ( q < 0.0 )
    {
      value = - value;
    }

  }

  return value;
}

double RelocSimulation::r8polyValueHorner(int m, double c[], double x){
  int i;
  double value;

  value = c[m];

  for ( i = m - 1; 0 <= i; i-- )
  {
    value = value * x + c[i];
  }

  return value;
}

double RelocSimulation::randomUniform(){
    //generate a random uniform number in the range [0,1)
    //double r = (double)rand() / RAND_MAX;
    double r = dis(rgen);
    return(r);
}

int RelocSimulation::getWardID(int ward){
    
    return((wards_pointer + ward)->wardnr);
}

double RelocSimulation::getWardArrivalRate(int ward){
    
    return((wards_pointer + ward)->arrivalRate);
}

double RelocSimulation::getWardServiceRate(int ward){
    
    return((wards_pointer + ward)->serviceRate);
}

vector<double> RelocSimulation::getWardRelocationProbabilities(int ward){
    
    return((wards_pointer + ward)->relocationProbabilities);
}

void RelocSimulation::calculateWardStateSpaceSize(int ward, int numberOfWards){
    
    return((wards_pointer + ward)->calculateWardStateSpace(numberOfWards));
}

int RelocSimulation::getWardStateSpaceSize(int ward){
    
    return((wards_pointer + ward)->wardStateSpaceSize);
}

int RelocSimulation::getWardCapacity(int ward){
    
    return((wards_pointer + ward)->capacity);
}