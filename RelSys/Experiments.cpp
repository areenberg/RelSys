
/* 
 * File:   Experiments.cpp
 * Author: anders
 * 
 * Created on January 26, 2021, 3:11 PM
 */

#include "Experiments.h"
#include "RelocSimulation.h"
#include "RelocEvaluation.h"
#include "QueueData.h"
#include "StatusBar.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;
using namespace std::chrono;


Experiments::Experiments(int nW, QueueData * wards, int sd):
wards_pointer(wards),
nWards(nW),
seed(sd)
{
    initialize();
}

Experiments::Experiments(const Experiments& orig) {
}

Experiments::~Experiments() {
}


void Experiments::initialize(){
    
    //set seed
    setSeed();
    
}



void Experiments::heuristicExperiment(int nRep, int widx, string label){
    
    cout << "Running experiment " << label << endl;
    cout << "Overall experiment seed: " << seed << endl; 
    
    //initialize vector of results
    distResults.resize(nRep);
    runTimeResults.resize(nRep,0);
    
    //seeds
    vector<int> seedList = generateSeeds(nRep);
    
    //EXPERIMENTS
    for (int rep=0; rep<seedList.size(); rep++){
        RelocEvaluation * relocEva = new RelocEvaluation[1];
        cout << "Rep " << (rep+1) << ", seed: " << seedList[rep] << endl;
        
        auto start = high_resolution_clock::now(); //start time
        
        //setup model object
        relocEva[0] = RelocEvaluation(nWards,wards_pointer);
        
        //simulate samples
        relocEva[0].runSimulation(seedList[rep],365,365,50);
    
        //evaluate ward
        //relocEva[0].simulateMarginalDist(50,10000); //simulate solution instead of solving numerically
//        vector<vector<int>> binMap = {{1,0,0,0,0}, //A
//                                      {0,1,0,0,0}, //B
//                                      {0,1,0,0,0}, //C
//                                      //{0,0,1,0,0}, //D
//                                      {0,1,0,0,0}, //E
//                                      {0,1,0,0,0}, //F
//                                      {0,1,0,0,0}, //G
//                                      {0,0,0,1,0}, //H
//                                      {0,0,0,1,0}, //I
//                                      {0,0,0,0,1}}; //J
//                                      //{0,0,1,0,0}};//K
//        relocEva[0].setBinMap(binMap);
        relocEva[0].runHeuristic(widx);
        
        //get runtime
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        runTimeResults[rep] = duration.count();
        if (rep==0){
            memResults = relocEva[0].vmemory; //estimated maximum virtual memory usage
        }
            
        //marginal distribution
        distResults[rep].resize(relocEva[0].marginalDist.size(),0);
        for (int j=0; j<relocEva[0].marginalDist.size(); j++){
            distResults[rep][j] = relocEva[0].marginalDist[j];
        }
        delete[] relocEva; //free model from memory
        
    }
    
    //save results
    saveResults(label);
}


void Experiments::simulationExperiment(double burnIn, double minTime,
        int nRep, int widx, string label, bool logNormal, double stdMult){
    
    cout << "Running experiment " << label << endl;
    cout << "Overall experiment seed: " << seed << endl; 
    
    //initialize vector of results
    distResults.resize(nRep);
    runTimeResults.resize(nRep,0);
    
    //seeds
    vector<int> seedList = generateSeeds(nRep);
    
    //EXPERIMENTS
    for (int rep=0; rep<seedList.size(); rep++){
        
        //setup model object
        RelocSimulation * sim_pointer = new RelocSimulation(nWards,wards_pointer);
        
        //run calculations
        sim_pointer->setSeed(seed);
        sim_pointer->disableTimeSampling();
        if (logNormal){
            sim_pointer->selectLogNormalServiceTime(stdMult);
        }
        auto start = high_resolution_clock::now();
        
        sim_pointer->simulate(burnIn,minTime);
        
        //get runtime
        auto stop = high_resolution_clock::now(); //stop time 
        auto duration = duration_cast<milliseconds>(stop - start); 
        runTimeResults[rep] = duration.count();
        if (rep==0){
            memResults = 0;
        }
        
        //marginal distribution
        distResults[rep].resize(sim_pointer->wardFreqDist[widx].size(),0);
        for (int j=0; j<sim_pointer->wardFreqDist[widx].size(); j++){
            distResults[rep][j] = sim_pointer->wardFreqDist[widx][j];
        }
        delete sim_pointer; //free model from memory
        cout << "rep. " << (rep+1) << " done." << endl;
    }
    cout << "All experiments done." << endl;
    //save results
    saveResults(label);    
}


void Experiments::KSExperiment(double burnIn, double minTime,
        int reps, bool logNormal, double stdMult, vector<vector<int>> obsDists){
    //Conducts a Kolmogorov-Smirnov bootstrap test on the observed
    //distributions in obsDists.
    //alpha is the significance level.
    //burnIn is the overall burn-in time.
    //minTime is the sim.time used in the reference distributions.
    //reps is the number of repetitions in the error distribution.
    //logNormal and stdMult defines the service time distribution.
    //obsDists are the observed (from collected data) marginal distributions.
    
    //(1) simulate reference distributions
    cout << "Starting KS test..." << endl;
    RelocSimulation* refSim = new RelocSimulation(nWards,wards_pointer);
    refSim->disableTimeSampling();    
    refSim->setSeed(seed);
    refSim->simulate(burnIn,minTime);
    
    vector<vector<double>> refDists(nWards);
    for (int widx=0; widx<nWards; widx++){
        refDists[widx].resize(refSim->wardFreqDist[widx].size(),0);
        for (int j=0; j<refSim->wardFreqDist[widx].size(); j++){
            refDists[widx][j] = (double)refSim->wardFreqDist[widx][j]/(double)refSim->nWardFreq[widx];
        }
        for (int j=1; j<refSim->wardFreqDist[widx].size(); j++){
            refDists[widx][j] += refDists[widx][j-1];
        }
    }
    delete refSim;
    
    //(2) simulate error samples
    cout << "Sampling error distribution..." << endl;
    vector<int> seedList = generateSeeds(reps);
    vector<int> maxWardSamples(nWards,0);
    vector<vector<double>> errorSamples(nWards);
    for (int widx=0; widx<nWards; widx++){
        for (int j=0; j<obsDists[widx].size(); j++){
            maxWardSamples[widx] += obsDists[widx][j];
        }
        errorSamples[widx].resize(reps,0);
    }
    
    StatusBar sbar(reps,30);
    for (int i=0; i<reps; i++){
        RelocSimulation* sim = new RelocSimulation(nWards,wards_pointer);
        sim->disableTimeSampling();
        sim->setSeed(seedList[i]);
        sim->simulate(burnIn,minTime,maxWardSamples);
        
        for (int widx=0; widx<nWards; widx++){
            vector<double> simDist(sim->wardFreqDist[widx].size(),0);
            for (int j=0; j<simDist.size(); j++){
                simDist[j] = (double)sim->wardFreqDist[widx][j]/(double)sim->nWardFreq[widx];
            }
            for (int j=1; j<simDist.size(); j++){
                simDist[j] += simDist[j-1];
            }
            errorSamples[widx][i] = KSStatistic(simDist,refDists[widx]);
        }
        
        if (i%10==0){
            double bval = i;
            sbar.updateBar(bval);
        }
        delete sim;
    }
    sbar.endBar();
    
    //(3) calculate test statistic and compare
    //to error distributions
    cout << "-----------------------" << endl;
    cout << "Results of KS test" << endl;
    cout << "-----------------------" << endl;         
    for (int widx=0; widx<nWards; widx++){
        cout << "WARD " << (widx+1) << ": ";
        vector<double> obsRelDist(obsDists[widx].size(),0);
        for (int j=0; j<obsRelDist.size(); j++){
            obsRelDist[j] = (double)obsDists[widx][j]/(double)maxWardSamples[widx];
        }
        for (int j=1; j<obsRelDist.size(); j++){
            obsRelDist[j] += obsRelDist[j-1];
        }
        cout << endl;
        double tst = KSStatistic(obsRelDist,refDists[widx]);
        double pval = (double)nLargerThan(tst,errorSamples[widx])/(double)errorSamples[widx].size();
        vector<double> confInt = wilsonScoreInterval(pval,reps);
        cout << "KS-test: " << tst << ", p-value: " << pval;
        if (pval>5e-2&&pval<1e-1){
            cout << " .";
        }else if(pval<5e-2){
            cout << " *";
            if (pval<1e-2){
                cout << "*";
            }
            if (pval<1e-3){
                cout << "*";
            }
        }
        cout << " (" << confInt[0] << "," << confInt[1] << ")";
        cout << endl;
    }
    
}

int Experiments::nLargerThan(double x, vector<double> v){
    //counts the number of elements in vector v that are larger
    //than x.
    int count=0;
    for (int i=0; i<v.size(); i++){
        if (v[i]>x){
            count++;
        }
    }
    return(count);
}

double Experiments::KSStatistic(vector<double>& dist1, vector<double>& dist2){
    //calculates the Kolmogorov-Smirnov test statistic using the
    //two CDFs dist1 and dist2.
    double df,mx = 0;
    for (int i=0; i<dist1.size(); i++){
        df = abs(dist1[i]-dist2[i]);
        if (df>mx){
            mx=df;
        }
    }
    return(mx);
}

void Experiments::chiSqExperiment(double burnIn, double minTime,
        int reps, bool logNormal, double stdMult,
        vector<vector<int>> obsDists, int method){
    //Conducts a Pearson's goodness of fit bootstrap test on the observed
    //distributions in obsDists.
    //method controls how the test statistic is calculated:
    //0: sum (o_ij-e_ij)^2/e_ij and 1: sum (o_ij-e_ij)^2 
    
    if (method==0){
        cout << "Starting method 0 test..." << endl;
        chiSqExperiment_method0(burnIn,minTime,
        reps,logNormal,stdMult,obsDists,method);
    }else if(method==1){
        cout << "Starting method 1 test..." << endl;
        chiSqExperiment_method1(burnIn,minTime,
        reps,logNormal,stdMult,obsDists,method);
    }else{
        cout << "Method type not recognized." << endl;
    }

}

void Experiments::chiSqExperiment_method0(double& burnIn, double& minTime,
        int& reps, bool& logNormal, double& stdMult,
        vector<vector<int>>& obsDists, int& method){
    
    //(1) simulate reference distributions and calculate expected frequencies
    
    //sum of samples in each ward
    vector<int> maxWardSamples(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        for (int j=0; j<obsDists[widx].size(); j++){
            maxWardSamples[widx] += obsDists[widx][j];
        }
    }
    
    RelocSimulation* refSim = new RelocSimulation(nWards,wards_pointer);
    refSim->disableTimeSampling();    
    refSim->setSeed(seed);
    refSim->simulate(burnIn,minTime);
    
    vector<vector<double>> eFreq_temp(nWards);
    for (int widx=0; widx<nWards; widx++){
        eFreq_temp[widx].resize(refSim->wardFreqDist[widx].size(),0);
        for (int j=0; j<refSim->wardFreqDist[widx].size(); j++){
            eFreq_temp[widx][j] = ((double)refSim->wardFreqDist[widx][j]/(double)refSim->nWardFreq[widx])*(double)maxWardSamples[widx];
        }
    }
    delete refSim;
    
    //consolidation indices in tails
    vector<vector<double>> eFreq(nWards);
    vector<vector<int>> consol_left(nWards);
    vector<vector<int>> consol_right(nWards);
    int j;
    double sm,lim=1;
    cout << "Bin consolidation:" << endl;
    for (int widx=0; widx<nWards; widx++){
        cout << "Ward " << (widx+1) << ": ";
        //left tail
        sm=0;
        j=0;
        cout << "Left: ";
        while (sm<=lim){
            cout << j << " ";
            consol_left[widx].push_back(j);
            sm += eFreq_temp[widx][j];
            j++;
        }
        //right tail
        sm=0;
        j=eFreq_temp[widx].size()-1;
        cout << "Right: ";
        while (sm<=lim){
            cout << j << " ";
            consol_right[widx].push_back(j);
            sm += eFreq_temp[widx][j];
            j--;
        }
        cout << endl;
        eFreq[widx] = consolidateDist(eFreq_temp[widx],
                consol_left[widx],consol_right[widx]);
        cout << "Dist: ";
        for (int l=0; l<eFreq[widx].size(); l++){
            cout << eFreq[widx][l] << " ";
        }
        cout << endl;
    }
    
    //(2) simulate error samples
    vector<vector<double>> errorSamples(nWards);
    for (int widx=0; widx<nWards; widx++){
        errorSamples[widx].resize(reps,0);
    }
    vector<int> seedList = generateSeeds(reps);
    StatusBar sbar(reps,30);
    double df;
    vector<double> simFreq,simFreq_temp;
    for (int i=0; i<reps; i++){
        RelocSimulation* sim = new RelocSimulation(nWards,wards_pointer);
        sim->disableTimeSampling();
        sim->setSeed(seedList[i]);
        sim->simulate(burnIn,minTime,maxWardSamples);
        
        for (int widx=0; widx<nWards; widx++){
            simFreq_temp.resize(sim->wardFreqDist[widx].size(),0);
            for (int j=0; j<simFreq_temp.size(); j++){
                simFreq_temp[j] = (double)sim->wardFreqDist[widx][j];
            }
            
            simFreq = consolidateDist(simFreq_temp,
                consol_left[widx],consol_right[widx]);
            for (int j=0; j<simFreq.size(); j++){
                df = simFreq[j]-eFreq[widx][j];
                errorSamples[widx][i] += pow(df,2.0)/eFreq[widx][j];
            }
        }
        
        if (i%10==0){
            double bval = i;
            sbar.updateBar(bval);
        }
        delete sim;
    }
    sbar.endBar();

    //(3) calculate test statistic and compare
    //to error distributions
    cout << "-----------------------" << endl;
    cout << "Results of test" << endl;
    cout << "-----------------------" << endl;
    double pval;
    vector<double> testStat(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        vector<double> obsFreq_temp(obsDists[widx].size(),0);
        for (int j=0; j<obsFreq_temp.size(); j++){
            obsFreq_temp[j] = (double)obsDists[widx][j];
        }
        vector<double> obsFreq = consolidateDist(obsFreq_temp,
            consol_left[widx],consol_right[widx]);

        for (int j=0; j<obsFreq.size(); j++){
            df = obsFreq[j]-eFreq[widx][j];
            testStat[widx] += pow(df,2.0)/eFreq[widx][j];
        }
        pval = (double)nLargerThan(testStat[widx],errorSamples[widx])/(double)errorSamples[widx].size();
        vector<double> confInt = wilsonScoreInterval(pval,reps);
        cout << "Ward " << (widx+1) << ": Test stat.: " << testStat[widx] << ", p-value: " << pval << " ("<< confInt[0] << "," << confInt[1] << ")" << endl;
    }
    double sumTest=0;
    for (int widx=0; widx<nWards; widx++){
        sumTest += testStat[widx];
    }
    vector<double> sumError(reps,0);
    for (int i=0; i<reps; i++){
        for (int widx=0; widx<nWards; widx++){
            sumError[i] += errorSamples[widx][i];
        }
    }
    pval = (double)nLargerThan(sumTest,sumError)/(double)sumError.size();
    vector<double> confInt = wilsonScoreInterval(pval,reps); 
    cout << "Overall test stat.: " << sumTest << ", p-value: " << pval << " ("<< confInt[0] << "," << confInt[1] << ")" << endl;
    
}

vector<double> Experiments::consolidateDist(vector<double>& dist,
        vector<int>& cleft, vector<int>& cright){
    
    int dim = dist.size()-((cleft.size()-1)+(cright.size()-1));
    vector<double> newDist(dim,0);
    
    //left and right consolidation
    for (int i=0; i<cleft.size(); i++){
        newDist[0] += dist[cleft[i]];
    }
    for (int i=0; i<cright.size(); i++){
        newDist[newDist.size()-1] += dist[cright[i]];
    }
    //fill in the rest
    int k = cleft.size();
    for (int i=1; i<(newDist.size()-1); i++){
        newDist[i] = dist[k]; 
        k++;
    }
    return(newDist);
}

void Experiments::chiSqExperiment_method1(double& burnIn, double& minTime,
        int& reps, bool& logNormal, double& stdMult,
        vector<vector<int>>& obsDists, int& method){

    //(1) simulate reference distributions and calculate expected frequencies
    
    //sum of samples in each ward
    vector<int> maxWardSamples(nWards,0);
    for (int widx=0; widx<nWards; widx++){
        for (int j=0; j<obsDists[widx].size(); j++){
            maxWardSamples[widx] += obsDists[widx][j];
        }
    }
    
    RelocSimulation* refSim = new RelocSimulation(nWards,wards_pointer);
    refSim->disableTimeSampling();    
    refSim->setSeed(seed);
    refSim->simulate(burnIn,minTime);
    
    vector<vector<double>> eFreq(nWards);
    for (int widx=0; widx<nWards; widx++){
        eFreq[widx].resize(refSim->wardFreqDist[widx].size(),0);
        for (int j=0; j<refSim->wardFreqDist[widx].size(); j++){
            eFreq[widx][j] = ((double)refSim->wardFreqDist[widx][j]/(double)refSim->nWardFreq[widx])*(double)maxWardSamples[widx];
        }
    }
    delete refSim;
    
    //(2) simulate error samples
    vector<vector<double>> errorSamples(nWards);
    for (int widx=0; widx<nWards; widx++){
        errorSamples[widx].resize(reps,0);
    }
    vector<int> seedList = generateSeeds(reps);
    StatusBar sbar(reps,30);
    double df;
    for (int i=0; i<reps; i++){
        RelocSimulation* sim = new RelocSimulation(nWards,wards_pointer);
        sim->disableTimeSampling();
        sim->setSeed(seedList[i]);
        sim->simulate(burnIn,minTime,maxWardSamples);
        
        for (int widx=0; widx<nWards; widx++){
            for (int j=0; j<sim->wardFreqDist[widx].size(); j++){
                df = (double)sim->wardFreqDist[widx][j]-eFreq[widx][j];
                errorSamples[widx][i] += pow(df,2.0);
            }
        }
        
        if (i%10==0){
            double bval = i;
            sbar.updateBar(bval);
        }
        delete sim;
    }
    sbar.endBar();
    
    //(3) calculate test statistic and compare
    //to error distributions
    cout << "-----------------------" << endl;
    cout << "Results of test" << endl;
    cout << "-----------------------" << endl;         
    vector<double> testStat(nWards,0);
    double pval;
    for (int widx=0; widx<nWards; widx++){
        
        for (int j=0; j<obsDists[widx].size(); j++){
            df = obsDists[widx][j]-eFreq[widx][j];
            testStat[widx] += pow(df,2.0);
        }
        
        pval = (double)nLargerThan(testStat[widx],errorSamples[widx])/(double)errorSamples[widx].size();
        vector<double> confInt = wilsonScoreInterval(pval,reps);
        cout << "Ward " << (widx+1) << ": Test stat.: " << testStat[widx] << ", p-value: " << pval << " ("<< confInt[0] << "," << confInt[1] << ")" << endl;
        
    }

    double sumTest=0;
    for (int widx=0; widx<nWards; widx++){
        sumTest += testStat[widx];
    }
    vector<double> sumError(reps,0);
    for (int i=0; i<reps; i++){
        for (int widx=0; widx<nWards; widx++){
            sumError[i] += errorSamples[widx][i];
        }
    }
    pval = (double)nLargerThan(sumTest,sumError)/(double)sumError.size();
    vector<double> confInt = wilsonScoreInterval(pval,reps); 
    cout << "Overall test stat.: " << sumTest << ", p-value: " << pval << " ("<< confInt[0] << "," << confInt[1] << ")" << endl;
    
}

vector<double> Experiments::wilsonScoreInterval(double p, int n){
    //computes binomial proportion confidence intervals
    //using the Wilson score interval method.
    
    double z = 1.959964;
    vector<double> intervals(2,0);
    
    double d0 = (1.0/(1.0+pow(z,2.0)/(double)n))*(p+(pow(z,2.0)/(2.0*(double)n)));
    double d1 = (z/(1.0+pow(z,2.0)/(double)n))*
    sqrt((p*(1-p))/n + (pow(z,4.0)/(4.0*pow((double)n,2.0))));
    
    //lower limit
    intervals[0] = d0-d1;
    //upper limit
    intervals[1] = d0+d1;
    
    return(intervals);
}


void Experiments::setSeed(){
    srand(seed);
}

double Experiments::randomUniform(){
    //generate a random uniform number in the range [0,1)
    return(rand()/((double) RAND_MAX));
}

int Experiments::randomInteger(int from, int to){
    int y = (int) round(randomUniform()*(to-from)+from); 
    return(y);
}

vector<int> Experiments::generateSeeds(int nRep){
    //generate nRep simulation seeds
    //for the experiments
    int y,from,to;
    vector<int> seedList(nRep,0);
    from = 100;
    to = 2147483647;
    
    for (int i=0; i<seedList.size(); i++){
        do {
            y = randomInteger(from,to);
        } while (numberExists(y,seedList));
        seedList[i] = y;
    }
    return(seedList);
}

bool Experiments::numberExists(int y, vector<int> list){
    //checks if the generated seed is already
    //on the list.
    
    bool exists = false;
    int i=0;
    while (i<list.size() && y!=list[i]){
        i++;
    }
    if (y==list[i]){
        exists = true;
    }
    
    return(exists);
}

void Experiments::saveResults(string label){
    //save experiment results in CSV-files.
    //distRes.csv contains the marginal distributions.
    //runMemRes.csv contains runtime and memory.
    
    //marginal distributions
    ofstream distResFile;
    string fName1 = "distRes_" + label + ".csv";
    distResFile.open(fName1);
    for(int i=0; i<distResults[0].size(); i++){
        for (int j=0; j<distResults.size(); j++){
            distResFile << distResults[j][i];
            if (j<(distResults.size()-1)){
                distResFile << "," << flush;
            }
        }
        distResFile << endl;
    }
    distResFile.close();
    
    //runtime and memory results
    ofstream runMemFile;
    string fName2 = "runMemRes_" + label + ".csv";
    runMemFile.open(fName2);
    for (int i=0; i<runTimeResults.size(); i++){
        runMemFile << runTimeResults[i];
        if (i<(runTimeResults.size()-1)){
            runMemFile << "," << flush;
        }
    }
    runMemFile << endl;
    runMemFile << memResults << endl;
    runMemFile.close();
    
    cout << "Results saved as " << fName1 << " and " << fName2 << endl;
    
}


