
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
checkAccuracy(false),
checkBurnIn(false),
burnInInit(true),
timeDepEnabled(false),        
stdMult(1.0),
accTol(5e-3),
accPref(true)        
{
    initializeSystem();
}

RelocSimulation::RelocSimulation(const RelocSimulation& orig) {
}

RelocSimulation::~RelocSimulation() {
}

void RelocSimulation::initializeSystem(){
    calculateArrivalRates();
    
    //pre-allocate memory for services and arrivals
    nextArrival = new Customer(0.0,0.0,0,0);
    maxOcc = 0;
    for (int widx=0; widx<nWards; widx++){
        maxOcc += getWardCapacity(widx);
    }
    service_array.resize(maxOcc);
    
    //set maximum generated number
    if (RAND_MAX==numeric_limits<int>::max()){
        mxRnd = RAND_MAX;
    }else{
        mxRnd = RAND_MAX+1;
    }
    
    //evaluate the (per bed) ward loads
    evalWardLoads();
    
}

void RelocSimulation::enableTimeDependency(QueuePerformance * qP){
    qPer=qP; //the queue occupancy distributions
    timeDepEnabled = true; //enable time-dependency
    cycleLen = wards_pointer->timeDep.size(); //get length of one cycle
}

void RelocSimulation::disableTimeSampling(){
    timeSamplingEnabled = false;
}

void RelocSimulation::setAccuracy(double a){
    //change the accuracy of the simulation.
    //only applicable when the simulation time is
    //set to -1
    accTol=a;
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

void RelocSimulation::evalWardLoads(){
    //evaluates bounds on the (per bed) ward loads
    
    double a;
    wardLoadUpperBounds.resize(nWards,0);
    wardLoadLowerBounds.resize(nWards,0);
    
    for (int widx=0; widx<nWards; widx++){
        wardLoadLowerBounds[widx]=(getWardArrivalRate(widx)/(getWardServiceRate(widx)*getWardCapacity(widx)));
        a=getWardArrivalRate(widx);
        for (int j=0; j<nWards; j++){
            if (widx!=j){
                a+=getWardArrivalRate(j)*getWardRelocationProbabilities(j)[widx];
            }
        }
        wardLoadUpperBounds[widx]=(a/(getWardServiceRate(widx)*getWardCapacity(widx)));
        //cout << "Ward " << (widx+1) << " load: [" << wardLoadLowerBounds[widx] << ";" << wardLoadUpperBounds[widx] << "]" << endl;
    }
    
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
    
    if (timeDepEnabled){
        updateArrivalTime_TimeDep(widx,pidx);
    }else{
        updateArrivalTime_Stationary(widx,pidx);
    }
    
}

void RelocSimulation::updateArrivalTime_TimeDep(int widx, int pidx){
    
    double delta=0.0;
    do{
        delta += randomExponential(arrivalRateMatrix[widx][pidx]);
    }while(randomUniform()>getWardTimeDep(widx,timeIndex(clock+delta)));
    
    nextArrivalTime[widx][pidx] += delta;
}

void RelocSimulation::updateArrivalTime_Stationary(int widx, int pidx){
    
    nextArrivalTime[widx][pidx] += randomExponential(arrivalRateMatrix[widx][pidx]);
    
}

int RelocSimulation::timeIndex(double cl){
    //return the time index corresponding to the
    //clock in cl
    return(floor(cl-floor(cl/cycleLen)*cycleLen));
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
        
    //create new patient arrival
    nextArrival->arrivalClock = nextArrivalTime[min_widx][min_pidx];
    nextArrival->serviceTime = genServiceTime(min_pidx);
    nextArrival->patientType = min_pidx;
    nextArrival->wardTarget = min_widx;
    nextArrival->active = true;
    nextArrival->serviceClock = numeric_limits<double>::max();
    
    updateArrivalTime(min_widx,min_pidx);
    
}


void RelocSimulation::simulate(double bIn, double minTime,
        vector<int> maxWardSamples, int minSamples){
    //burnIn is the burn-in time.
    //minTime is the minimal value of the clock.
    //minSamples is the minimal number of sampled open/blocked times.
    //default value for minSamples=50.
    
    //variables
    int targetWard;
    bool succeeded;

    //--------------------
    //Initialize
    //--------------------
    
    //for evaluating time-out
    timeOut=false;
    
    capUse.resize(nWards,0);
    wardOccupancy.resize(nWards);
    nextArrivalTime.resize(nWards);
    wilsonIntervals.resize(2,0);
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
    if (bIn>=0){
        burnIn = bIn;
    }else{
        burnIn = numeric_limits<double>::max();
        bInSize=30;
        burnInSamples.resize(2);
        burnInSamples[0].resize(bInSize);
        burnInSamples[1].resize(bInSize);
        bInOrder = {0,1};
        clockDis=0.0;
        disIdx=0;
        checkBurnIn = true;
    }
    //simulation time
    if (minTime>=0){
        simTime = minTime;
    }else{
        simTime=numeric_limits<double>::max();
        checkAccuracy=true;
    }
    
    //simulation clock
    clock = 0;
    
    //patient arrivals
    initializeArrivalTimes(clock);
    generateArrival();
    
    //patients in service
    inService = 0;
    serIdx = 0;
    
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
            maxWrdSam[i] = numeric_limits<int>::max()-1;
        }
    }else{
        for (int i=0; i<nWards; i++){
            maxWrdSam[i] = maxWardSamples[i];
        }
    }
    //minSamples cannot be less than 2
    if (timeSamplingEnabled && minSamples<2){
        minSamples = 2;
    }
    
    //-------------------
    //Simulate
    //-------------------
    
    auto start = chrono::system_clock::now();
    cout << "Running simulation..." << flush;
    while (timeOut==false && ((accuracy()>accTol && clock<simTime && wardSamplesToGo()>0) || (timeSamplingEnabled && minTimeSamples()<minSamples)) ){
        
        
        if (inService>0){
            
            nextServiceIdx();
            
            if (service_array[serIdx].serviceClock<nextArrival->arrivalClock){
            
                clock = service_array[serIdx].serviceClock;
                targetWard = service_array[serIdx].wardTarget;
                attemptDischarge();
                if (timeSamplingEnabled && succeeded){
                    blockedTimeTracking(targetWard);
                }
                if (checkBurnIn){
                    evaluateBurnIn();
                }
                
            }else{

                clock = nextArrival->arrivalClock;
                attemptAdmission(succeeded);
                if (timeSamplingEnabled && succeeded){
                    openTimeTracking(nextArrival->wardTarget);
                }
                generateArrival();
            }

        }else{
            clock = nextArrival->arrivalClock; 
            attemptAdmission(succeeded);
            if (timeSamplingEnabled && succeeded){
                openTimeTracking(nextArrival->wardTarget);
            }
            generateArrival();

        }

        interElapsed = chrono::duration_cast<chrono::seconds>(chrono::system_clock::now() - start).count();
        checkTerminate();
    }
    cout << "done." << endl;

    performanceMeasures();
    if (timeSamplingEnabled){
        //re-sample time tracking
        subsetTimeSamples(minSamples);
        //print open/blocked time sample sizes
        //printTimeSamples();
    }
    
    auto stop = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Runtime of simulation: " << elapsed.count() << " milliseconds\n" << endl;
}

void RelocSimulation::evaluateBurnIn(){
    
    double tDiff = clock-clockDis;
    clockDis = clock;
    
    if (burnInInit && disIdx<bInSize){
        burnInSamples[0][disIdx]=tDiff;
        disIdx++;
    }else if(burnInInit){
        burnInInit=false;
        disIdx=0;
    }else{
        burnInSamples[bInOrder[1]][disIdx]=tDiff;
        disIdx++;
    }
    
    //evaluate
    if (disIdx==bInSize && burnInInit==false){
        
        bool larger=wilcoxonRankSum(burnInSamples[bInOrder[0]],burnInSamples[bInOrder[1]]);
        
        if (larger){
            disIdx=0;
            int xx=bInOrder[0];
            bInOrder[0]=bInOrder[1];
            bInOrder[1]=xx;
        }else{
            checkBurnIn=false;
            burnIn=clock;
            cout << " simulation burn-in at: " << burnIn << " (continuing simulation) ... " << flush;
        }
    }
    
}

void RelocSimulation::checkTerminate(){
    //decides if the simulation should be terminated
    
    if (interElapsed>1 && (timeSamplingEnabled||checkAccuracy)){ //check after 1 second
        
        for (int widx=0; widx<nWards; widx++){
            if (timeSamplingEnabled && openTimes[widx].empty() && //terminate early if no hope to obtain a sufficient number of open/shortage time-samples
                    wardLoadLowerBounds[widx]<0.1){
                cout << "Queue " << (widx+1) << " load caused simulation to terminate early." << flush;
                timeOut=true;
            }else if(checkAccuracy && interElapsed>200 && ((!accPref && nWardFreq[widx]<100) || (accPref && nWardFreqPref[widx]<100))){
                skipAccuracy.push_back(widx);
            }
        }
        
    }
    
}

int RelocSimulation::wardSamplesToGo(){
    //maximum number of samples in the marginal
    //frequency distributions that haven't been
    //collected yet.
    int df, mx=0;
    for (int widx=0; widx<nWards; widx++){
        if (maxWrdSam[widx]==numeric_limits<int>::max()){
            return(numeric_limits<int>::max());
        }
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
    wardDenDist.resize(nWards);//marginal density distribution
    nWardFreq.resize(nWards); //samples in marginal freq. dist.
    wardFreqDistPref.resize(nWards);//marginal frequency distribution
    wardDenDistPref.resize(nWards);//marginal density distribution
    nWardFreqPref.resize(nWards); //samples in marginal freq. dist.
    freqDist.resize(nWards); //frequency distribution
    
    for (int widx=0; widx<nWards; widx++){
        int c = getWardCapacity(widx)+1;
        wardFreqDist[widx].resize(c,0);
        wardDenDist[widx].resize(c,0);
        nWardFreq[widx] = 0;
        wardFreqDistPref[widx].resize(c,0);
        wardDenDistPref[widx].resize(c,0);
        nWardFreqPref[widx] = 0;
        freqDist[widx].resize(nWards);
        for (int pidx=0; pidx<nWards; pidx++){
            freqDist[widx][pidx].resize(c,0);
        }
    }
    
}

void RelocSimulation::occupancyDistTracking(int &targetWard, int &patientType){
    
    if (clock>burnIn){
        
        //distribution observed only for preferred arrivals
        //note: preferred arrivals observe all wards at the same time
        if (targetWard==patientType){
            for (int widx=0; widx<nWards; widx++){
                if (nWardFreqPref[widx]<maxWrdSam[widx]){
                    wardFreqDistPref[widx][capUse[widx]]++;
                    nWardFreqPref[widx]++;
                }
            }
        }
        
        //distribution observed from all arrivals to the target ward (both preferred and alternative)
        if ( (targetWard!=patientType && getWardCapacity(patientType)==capUse[patientType]) || 
                targetWard==patientType ){
            //record distribution over all patient types in the ward
            for (int pidx=0; pidx<nWards; pidx++){
                freqDist[targetWard][pidx][wardOccupancy[targetWard][pidx]]++;
            }
            
            if (nWardFreq[targetWard]<maxWrdSam[targetWard]){ //check sampling does not exceed maximum sampling limit
                //wardFreqDist[targetWard][wardSum]++;
                wardFreqDist[targetWard][capUse[targetWard]]++;
                nWardFreq[targetWard]++;
            }
            
        }
    }
}

double RelocSimulation::accuracy(){
    //calculates the accuracy of the density
    //distributions using Wilson score intervals
    //note: high accuracy is indicated with a small value
    
    if (clock<burnIn||checkAccuracy==false){
        return(numeric_limits<double>::max());
    }else{
        wilsonMax = -1;
        for (int widx=0; widx<nWards; widx++){
            
            if (!skipWardAccuracy(widx)){
            
                if ((!accPref && nWardFreq[widx]==0) || (accPref && nWardFreqPref[widx]==0)){
                    return(numeric_limits<double>::max());
                }
                
                for (int j=0; j<wardFreqDist[widx].size(); j++){
                
                    if (!accPref){
                        wardDenDist[widx][j] = (double)wardFreqDist[widx][j]/(double)nWardFreq[widx];
                    }else{
                        wardDenDistPref[widx][j] = (double)wardFreqDistPref[widx][j]/(double)nWardFreqPref[widx];
                    }

                    wilsonScoreInterval(wilsonSpan,j,widx);
                    if (wilsonSpan>wilsonMax){
                        wilsonMax=wilsonSpan;
                    } 
                }
            
            }
            
        }
        return(wilsonMax);
    }
}

void RelocSimulation::setAccPref(bool ap){
    //set the sampling method in auto simulation time
    //true: all preferred arrivals, false: all arrivals to each ward
    accPref = ap;
}

bool RelocSimulation::skipWardAccuracy(int &widx){
    
    if (skipAccuracy.empty()){
        return(false);
    }else{
    
        int i=0;
        while(i<skipAccuracy.size() && widx!=skipAccuracy[i]){
            i++;
        }
        if (i==skipAccuracy.size()){
            return(false);
        }else{
            return(true);
        }
    }
    
}

void RelocSimulation::performanceMeasures(){
    //derive the fundamental performance measures
    
    //derive estimate of density and expected occupancy
    wardDenDist.resize(nWards);
    wardDenDistPref.resize(nWards);
    expectedOccupancy.resize(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        wardDenDist[widx].resize(wardFreqDist[widx].size(),0);
        wardDenDistPref[widx].resize(wardFreqDistPref[widx].size(),0);

        //density distribution for preferred arrivals
        if (nWardFreqPref[widx]==0){
           wardDenDistPref[widx][0]=1.0;
           expectedOccupancy[widx]=0.0;
        }else{
            for (int i=0; i<wardFreqDistPref[widx].size(); i++){
                wardDenDistPref[widx][i] = (double)wardFreqDistPref[widx][i]/(double)nWardFreqPref[widx];
                expectedOccupancy[widx] += i*wardDenDistPref[widx][i];
            }
        }

        //density distribution for all arrivals to the ward
        if (nWardFreq[widx]==0){
           wardDenDist[widx][0]=1.0;
        }else{
            for (int i=0; i<wardFreqDist[widx].size(); i++){
                wardDenDist[widx][i] = (double)wardFreqDist[widx][i]/(double)nWardFreq[widx];
            }
        }
        
    }
    //blocking probability and fraction of occupied servers
    blockingProbability.resize(nWards,0);
    blockingProbabilityPref.resize(nWards,0);
    expOccFraction.resize(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        blockingProbability[widx] = wardDenDist[widx][wardDenDist[widx].size()-1];
        blockingProbabilityPref[widx] = wardDenDistPref[widx][wardDenDistPref[widx].size()-1];
        expOccFraction[widx] = expectedOccupancy[widx]/getWardCapacity(widx);
    }
    
}


void RelocSimulation::openTimeTracking(int &targetWard){
    
    if (capUse[targetWard]==getWardCapacity(targetWard)){
        if (clock>burnIn){
            double sample = clock - wardStateClocks[targetWard];
            openTimes[targetWard].push_back(sample);
            nOpenTimeSamples[targetWard]++;
        }
        wardStateClocks[targetWard] = clock;
    }
    
}
    
void RelocSimulation::blockedTimeTracking(int &targetWard){
    
    if (capUse[targetWard]== (getWardCapacity(targetWard)-1) ){
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
    
    //adjust use of capacity in ward
    capUse[service_array[serIdx].wardTarget]--;
    wardOccupancy[service_array[serIdx].wardTarget][service_array[serIdx].patientType]--;
    
    //track occupancy if time dep. enabled
    if (timeDepEnabled){
        if (clock>burnIn){
            qPer->discharge(clock,service_array[serIdx].wardTarget,true);
        }else{
            qPer->discharge(clock,service_array[serIdx].wardTarget,false);
        }
            
    }
    
    return(true);
}

void RelocSimulation::attemptAdmission(bool &succeeded){

    //track regardless of admission accepted or rejected
    occupancyDistTracking(nextArrival->wardTarget,
                    nextArrival->patientType);
    
    if ( (nextArrival->wardTarget==nextArrival->patientType &&
            capUse[nextArrival->wardTarget]<getWardCapacity(nextArrival->wardTarget)) || (nextArrival->wardTarget!=nextArrival->patientType &&
            capUse[nextArrival->patientType]==getWardCapacity(nextArrival->patientType) &&
            capUse[nextArrival->wardTarget]<getWardCapacity(nextArrival->wardTarget)) ){
        
        //find an available slot
        insIdx=0;
        while (service_array[insIdx].active){
            insIdx++;
        }

        //insert into service array
        service_array[insIdx] = Customer(nextArrival->arrivalClock,nextArrival->serviceTime,
                nextArrival->wardTarget,nextArrival->patientType);
        //calculate and insert clock at discharge
        service_array[insIdx].serviceClock = clock + service_array[insIdx].serviceTime;
        
        //adjust use of capacity in the ward
        capUse[service_array[insIdx].wardTarget]++;
        wardOccupancy[service_array[insIdx].wardTarget][service_array[insIdx].patientType]++;
        
//        cout << "admission: " << insIdx << endl;
                
        //adjust number of patients currently in service
        inService++;
        
        //track occupancy if time dep. enabled
        if (timeDepEnabled){
            if (clock>burnIn){
                qPer->arrival(clock,service_array[insIdx].wardTarget,true);
            }else{
                qPer->arrival(clock,service_array[insIdx].wardTarget,false);
            }
            
        }
        
        succeeded = true;
    }else{
        succeeded = false;
    }
    
//    return(Ok);
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

void RelocSimulation::wilsonScoreInterval(double &wilsonSpan, int &j, int &widx){
    //calculates binomial proportion confidence intervals
    //using the Wilson score interval method.
    
    if (!accPref){
        wilsonSpan = 2*(( 1.959964/( 1.0+pow(1.959964,2.0) /(double)nWardFreq[widx]))*
        sqrt((wardDenDist[widx][j]*(1-wardDenDist[widx][j]))/nWardFreq[widx] + (pow(1.959964,4.0)/(4.0*pow((double)nWardFreq[widx],2.0)))));
    }else{
        wilsonSpan = 2*(( 1.959964/( 1.0+pow(1.959964,2.0) /(double)nWardFreqPref[widx]))*
        sqrt((wardDenDistPref[widx][j]*(1-wardDenDistPref[widx][j]))/nWardFreqPref[widx] + (pow(1.959964,4.0)/(4.0*pow((double)nWardFreqPref[widx],2.0)))));
    }

}

bool RelocSimulation::wilcoxonRankSum(vector<double> x, vector<double> y){
    //conducts the Wilcoxon rank-sum test.
    //returns true if samples x are stochastically larger
    //than samples y - i.e. the null hypothesis that the samples are from the
    //same distribution is rejected.
    
//    double mn=0;
//    for (int i=0; i<x.size(); i++){
//        mn += x[i];
//    }
//    mn /= x.size();
//    cout << "x mean: " << mn << endl;
//    mn=0;
//    for (int i=0; i<y.size(); i++){
//        mn += y[i];
//    }
//    mn /= y.size();
//    cout << "y mean: " << mn << endl;
    
    int u=0;
    for (int i=0; i<x.size(); i++){
        for (int j=0; j<y.size(); j++){
            if (x[i]>y[j]){
                u++;
            }
        }
    }
    
    double m_u = (x.size()*y.size())/2.0;
    double s_u = sqrt((x.size()*y.size()*(x.size()+y.size()+1.0))/12.0);
    double z = (u-m_u)/s_u;
//    cout << "u=" << u << " m_u=" << m_u << " s_u=" << s_u << " z=" << z << endl;
    double zCrit = 1.036433; //critical z-score value (one-tailed test with alpha=0.15)
    if (z>zCrit){
        return(true); //null hypothesis rejected
    }else{
        return(false); //null hypothesis accepted (the populations are equal)
    }
    
}

double RelocSimulation::randomUniform(){
    //generate a random uniform number in the range [0,1)
    return((double)rand()/((double)mxRnd));
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

double RelocSimulation::getWardTimeDep(int ward, int timeIndex){
    
    return((wards_pointer + ward)->timeDep[timeIndex]);
}