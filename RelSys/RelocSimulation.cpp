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
#include "Customer.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <limits.h>
#include <cstdlib>
#include <chrono>

using namespace std;

RelocSimulation::RelocSimulation(int nW, QueueData * wards):
wards_pointer(wards),
nWards(nW),
serTimeExponential(true),
timeSamplingEnabled(true),
stdMult(1.0)        
{
    initializeSystem();
}

RelocSimulation::RelocSimulation(const RelocSimulation& orig) {
}

RelocSimulation::~RelocSimulation() {
}

void RelocSimulation::initializeSystem(){
    calculateArrivalRates();
    
    //pre-allocate memory for services
    nextArrival.resize(1);
    maxOcc = 0;
    for (int widx=0; widx<nWards; widx++){
        maxOcc += getWardCapacity(widx);
    }
    service_array.resize(maxOcc);
}

void RelocSimulation::disableTimeSampling(){
    timeSamplingEnabled = false;
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
    simSeed = seed;
    srand(simSeed);
    //mt19937 rgen(simSeed);
    //default_random_engine rgen{static_cast<long unsigned int>(seed)};
}

void RelocSimulation::initializeArrivalTimes(double currentClock){
    //allocate memory and fill in the entire nextArrivalTime matrix
    for (int i=0; i<nWards; i++){
        for (int j=0; j<nWards; j++){
            nextArrivalTime[i][j]=currentClock;
        }
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

void RelocSimulation::generateArrival(){
    
    //find minimum arrival time
    double mn = numeric_limits<double>::max(); 
    for (int widx=0; widx<nWards; widx++){
       for (int pidx=0; pidx<nWards; pidx++){
           if (nextArrivalTime[widx][pidx]<mn){
               mn = nextArrivalTime[widx][pidx];
               min_widx = widx; min_pidx = pidx;
            }
        }
    }        
        
    //create
    nextArrival[0] = Customer(nextArrivalTime[min_widx][min_pidx],
            genServiceTime(min_pidx),min_widx,min_pidx);
    updateArrivalTime(min_widx,min_pidx);
    
}

//void RelocSimulation::generateArrivalList(double currentClock){
//    //generate a long list of future arrivals with arrival times offset
//    //by currentClock
//    
//    initializeArrivalTimes(currentClock);
//    
//    //append patients to the list
//    double arrClock, serTime, mn;
//    int min_widx, min_pidx;
//    //vector<int> genWidx(nWards,0);
//    //vector<int> genPidx(nWards,0);
//    
//    for (int i=0; i<patientArraySize; i++){
//        
//        mn = numeric_limits<double>::max();
//        
//        //find minimum arrival time        
//        for (int widx=0; widx<nWards; widx++){
//            for (int pidx=0; pidx<nWards; pidx++){
//                if (nextArrivalTime[widx][pidx]<mn){
//                    mn = nextArrivalTime[widx][pidx];
//                    min_widx = widx; min_pidx = pidx;
//                }
//            }
//        }        
//        
//        
//        //append patient
//        arrClock = nextArrivalTime[min_widx][min_pidx];
//        serTime = genServiceTime(min_pidx);
//        arrival_array[i] = Customer(arrClock,serTime,min_widx,min_pidx);
//        
//        //update ward-patient
//        updateArrivalTime(min_widx,min_pidx);
//    }
//    
//}

void RelocSimulation::simulate(double bIn, double minTime,
        vector<int> maxWardSamples, int minSamples){
    //burnIn is the burn-in time.
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
//    for (int i=0; i<nWards; i++){
//        cout << getWardCapacity(i) << " " << flush;
//    }
//    cout << endl;

    
    //variables
    int arrIdx, targetWard;
    bool succeeded;

    //--------------------
    //Initialize
    //--------------------
    
    capUse.resize(nWards,0);
    wardOccupancy.resize(nWards);
    nextArrivalTime.resize(nWards);
    for (int widx=0; widx<nWards; widx++){
         wardOccupancy[widx].resize(nWards,0);
         nextArrivalTime[widx].resize(nWards,0);
    }
    
    nOpenTimeSamples.resize(nWards,0);
    nBlockedTimeSamples.resize(nWards,0);
    wardStateClocks.resize(nWards,0);
    
    //dummy customers in service_array
    for (int i=0; i<maxOcc; i++){
        service_array[i] = Customer(numeric_limits<double>::max(),0.0,-1,-1);
        service_array[i].active=false;
    }
    
    //burn-in time
    burnIn = bIn;
    
    //simulation clock
    clock = 0;
    
    //patient arrivals
    //generateArrivalList(clock);
    //arrIdx = 0;
    initializeArrivalTimes(clock);
    generateArrival();
    
    //patients in service
    inService = 0;
    serIdx = 0;
    
    //occupancy of the system
    vector<int> capUse(nWards,0);
    
    //tracking
    openTimes.resize(nWards);
    blockedTimes.resize(nWards);
    //initialize the frequency/density
    //ward-patient distributions.
    initFreqDenDist();  
    //configure ward sample limit
    maxWrdSam.resize(nWards);
    if (maxWardSamples.size()==1 && maxWardSamples[0]==-1){
        for (int i=0; i<nWards; i++){
            maxWrdSam[i] = numeric_limits<int>::max();
        }
    }else{
        for (int i=0; i<nWards; i++){
            maxWrdSam[i] = maxWardSamples[i];
        }
    }
    
    //-------------------
    //Simulate
    //-------------------
    auto start = chrono::system_clock::now();
    cout << "Running simulation..." << flush;
    while ( (clock<minTime && wardSamplesToGo()>0) || (timeSamplingEnabled && minTimeSamples()<minSamples) ){
        
        if (inService>0){
            
            nextServiceIdx();
            
            if (service_array[serIdx].serviceClock<nextArrival[0].arrivalClock){
            
                clock = service_array[serIdx].serviceClock;
                targetWard = service_array[serIdx].wardTarget;
                succeeded = attemptDischarge();
                if (timeSamplingEnabled && succeeded){
                    blockedTimeTracking(targetWard);
                }
                
            }else{

                clock = nextArrival[0].arrivalClock;
                succeeded = attemptAdmission(arrIdx);
                if (timeSamplingEnabled && succeeded){
                    openTimeTracking(nextArrival[0].wardTarget);
                }
                generateArrival();
            }

        }else{
            clock = nextArrival[0].arrivalClock; 
            succeeded = attemptAdmission(arrIdx);
            if (timeSamplingEnabled && succeeded){
                openTimeTracking(nextArrival[0].wardTarget);
            }
            generateArrival();

        }
        
        updateOccupancy();
        
        //generate a new batch of patient arrivals
        //if (arrIdx==(patientArraySize-1)){
            //generateArrivalList(clock);
            //arrIdx = 0;
        //}
        
    }
    cout << "done." << endl;
    freqToDensity();
    performanceMeasures();
    if (timeSamplingEnabled){
        //re-sample time tracking
        subsetTimeSamples(minSamples);
        //print open/blocked time sample sizes
        printTimeSamples();
    }
    
    auto stop = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(stop - start);
    runtime = elapsed.count();
    cout << "runtime: " << runtime << " ms" << endl;
}

int RelocSimulation::wardSamplesToGo(){
    //maximum number of samples in the marginal
    //frequency distributions that haven't been
    //collected yet.
    int df, mx=0;
    for (int widx=0; widx<nWards; widx++){
        df = maxWrdSam[widx]-nWardFreq[widx];
        if (df>mx){
            mx=df;
        }
    }
    return(mx);
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

double RelocSimulation::genServiceTime(int idx){
    
    double y;
    
    if (serTimeExponential){
        //generate random exponential double
        y = randomExponential(getWardServiceRate(idx));
    }else{
        //generate random log-normal double
        double mean = 1.0/getWardServiceRate(idx);
        double std = stdMult/getWardServiceRate(idx);
        y = randomLogNormal(mean,std);
    }
    
    return(y);


}

void RelocSimulation::selectLogNormalServiceTime(double mult){
    
    serTimeExponential = false;
    stdMult = mult;
    cout << "Service time switched to log-normal with multiplier " << stdMult << endl;
    
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
    cout << "*Included and total open-time samples*" << endl;
    for (int i=0; i<openTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << openTimes[i].size() << " (" << nOpenTimeSamples[i] << ")";
        if (i<(openTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << "\n" << endl;
    cout << "*Included and total blocked-time samples*" << endl;
    for (int i=0; i<blockedTimes.size(); i++){
        cout << "Ward" << (i+1) << ": " << blockedTimes[i].size() << " (" << nBlockedTimeSamples[i] << ")";
        if (i<(blockedTimes.size()-1)){
            cout << ", ";
        }
    }
    cout << endl;
    cout << "----------------------------\n" << endl;
    
}

void RelocSimulation::initFreqDenDist(){
    
    wardFreqDist.resize(nWards);//marginal frequency distribution
    nWardFreq.resize(nWards); //samples in marginal freq. dist.
    freqDist.resize(nWards); //frequency distribution
    denDist.resize(nWards); //density distribution
    
    for (int widx=0; widx<nWards; widx++){
        int c = getWardCapacity(widx)+1;
        wardFreqDist[widx].resize(c,0);
        nWardFreq[widx] = 0;
        freqDist[widx].resize(nWards);
        denDist[widx].resize(nWards);
        for (int pidx=0; pidx<nWards; pidx++){
            freqDist[widx][pidx].resize(c,0);
            denDist[widx][pidx].resize(c,0);
        }
    }
    
}

void RelocSimulation::updateOccupancy(){
    
    //ward occupancy (all patients)
    for (int i=0; i<nWards; i++){
        capUse[i]=0;
    }
    
    for (int i=0; i<maxOcc; i++){
        if (service_array[i].active){
            capUse[service_array[i].wardTarget]++;
        }
    }
    
//    for (int i=0; i<nWards; i++){
//        cout << capUse[i] << " " << flush;
//    }
//    cout << endl;
    
    //ward occupancy (each patient type)
    for (int widx=0; widx<nWards; widx++){
        for (int i=0; i<nWards; i++){
            wardOccupancy[widx][i]=0;
        }
    }
    for (int i=0; i<maxOcc; i++){
        if (service_array[i].active){
            wardOccupancy[service_array[i].wardTarget][service_array[i].patientType]++;
        }
    }
    
}

void RelocSimulation::occupancyDistTracking(int &targetWard, int &patientType){
    
    if (clock>burnIn){
        if ( (targetWard!=patientType && getWardCapacity(patientType)==capUse[patientType]) ||
                targetWard==patientType ){
            //record distribution over all patient types in the ward
            int wardSum=0;
            for (int pidx=0; pidx<nWards; pidx++){
                freqDist[targetWard][pidx][wardOccupancy[targetWard][pidx]]++;
                wardSum += wardOccupancy[targetWard][pidx];
            }
            
            if (nWardFreq[targetWard]<maxWrdSam[targetWard]){ //check sampling does not exceed maximum sampling limit
                wardFreqDist[targetWard][wardSum]++;
                nWardFreq[targetWard]++;
            }
            
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

void RelocSimulation::performanceMeasures(){
    //derive the fundamental performance measures
    
    //derive estimate of density and expected occupancy
    wardDenDist.resize(nWards);
    expectedOccupancy.resize(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        wardDenDist[widx].resize(wardFreqDist[widx].size(),0);
        if (nWardFreq[widx]==0){
           wardDenDist[widx][0]=1.0;
           expectedOccupancy[widx]=0.0;
        }else{
            for (int i=0; i<wardFreqDist[widx].size(); i++){
                wardDenDist[widx][i] = (double)wardFreqDist[widx][i]/(double)nWardFreq[widx];
                expectedOccupancy[widx] += i*wardDenDist[widx][i];
            }
        }
    }
    //blocking probability and fraction of occupied servers
    blockingProbability.resize(nWards,0);
    expOccFraction.resize(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        blockingProbability[widx] = wardDenDist[widx][wardDenDist[widx].size()-1];
        expOccFraction[widx] = expectedOccupancy[widx]/getWardCapacity(widx);
    }
    
}


void RelocSimulation::openTimeTracking(int &targetWard){
    
    if ((capUse[targetWard]+1)==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
            openTimes[targetWard].push_back(sample);
            nOpenTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}
    
void RelocSimulation::blockedTimeTracking(int &targetWard){
    
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



bool RelocSimulation::attemptDischarge(){
    
    //subtract patient from service
    inService--;
    
    service_array[serIdx].active=false;
//    cout << "served: " << serIdx << endl;
    
    return(true);
}

bool RelocSimulation::attemptAdmission(int &arrIdx){
    bool Ok = false;
    if (nextArrival[0].wardTarget==nextArrival[0].patientType &&
            capUse[nextArrival[0].wardTarget]<getWardCapacity(nextArrival[0].wardTarget)){
        Ok = true;
    }else if(nextArrival[0].wardTarget!=nextArrival[0].patientType &&
            capUse[nextArrival[0].patientType]==getWardCapacity(nextArrival[0].patientType) &&
            capUse[nextArrival[0].wardTarget]<getWardCapacity(nextArrival[0].wardTarget)){
        Ok = true;
    }
    
    if (Ok){
        
        //find an available slot
        insIdx=0;
        while (service_array[insIdx].active){
            insIdx++;
        }

        //insert into service array
        service_array[insIdx] = Customer(nextArrival[0].arrivalClock,nextArrival[0].serviceTime,
                nextArrival[0].wardTarget,nextArrival[0].patientType);
        //calculate and insert clock at discharge
        service_array[insIdx].serviceClock = clock + service_array[insIdx].serviceTime;
        
//        cout << "admission: " << insIdx << endl;
                
                
        //adjust number of patients currently in service
        inService++;
    }
    
    //track regardless of admission accepted or rejected
    occupancyDistTracking(nextArrival[0].wardTarget,
                    nextArrival[0].patientType);
    
    //move to next arrival
    //arrIdx++;
    
    return(Ok);
}

void RelocSimulation::nextServiceIdx(){
    
    double mn = numeric_limits<double>::max();
    for (int i=0; i<maxOcc; i++){
        if (service_array[i].active && service_array[i].serviceClock<mn){
            mn = service_array[i].serviceClock;
            serIdx = i;
         }
    }
}

int RelocSimulation::minTimeSamples(){
        int mn = numeric_limits<int>::max();
        for (int i=0; i<openTimes.size(); i++){
            if (openTimes[i].size()<mn){
                mn = openTimes[i].size();
            }
        }
        for (int i=0; i<blockedTimes.size(); i++){
            if (blockedTimes[i].size()<mn){
                mn = blockedTimes[i].size();
            }
        }
        //cout << "done." << endl;
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
    return(rand()/((double) RAND_MAX));
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