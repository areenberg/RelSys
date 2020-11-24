//MODIFIED VERSION OF EMPHT
//Empht is a C-program developed by: S. Asmussen, O. Nerman & M. Olsson
//Link to the original Empht source code: 
//https://web.archive.org/web/20180617130551/http://home.math.au.dk/asmus/pspapers.html


/* 
 * File:   PhaseFitter.h
 * Author: Anders Reenberg Andersen
 *
 * Created on September 1, 2020, 3:26 PM
 */

#ifndef PHASEFITTER_H
#define PHASEFITTER_H


#include <vector>

using namespace std;

class PhaseFitter {
public:
    PhaseFitter();
    PhaseFitter(const PhaseFitter& orig);
    virtual ~PhaseFitter();
    
    //methods
    void run(int NoOfEMsteps=100, int intSeed=123, int inDType=6, double trP=100, double dt=1e-3);
    
    //type of output
    void setGeneralPH(int p);
    void setHyberExponential(int p);
    void setSumOfExponentials(int p);
    void setCoxian(int p);
    void setCoxianGeneral(int p);
    
    //type of input
    void setInputSample(vector<double> &samples);
    void setInputDensity(vector<vector<double>> &q_in, vector<double> &p_in);

    //the final fitted parameters
    vector<double> init_dist;
    vector<double> exit_rate_vector;
    vector<vector<double>> ph_generator; 
    
private:
    
    //variables and pointers
    double *obs, *weight, *censur, *cweight, *lower, *upper, *intweight;
    double  SumOfWeights, SumOfCensored, SumOfInt;
    int  *partner, NoOfObs, NoOfCensored, NoOfInt; 
    
    int phases, sampleOrDensity, densityType, PHDistribution, integerSeed;
    double truncationPoint, densityInterval;   
    
    double * pi_in_point;
    
    vector<double> sample_vector; //contains all samples
    vector<double> weight_vector; //contains all sample weights (usually 1)
    
    //private EMpht methods
    double *v_alloc(int n);
    int *int_v_alloc(int n);
    double **m_alloc(int n, int k);
    int **int_m_alloc(int n, int k);
    double ***m3_alloc(int n, int k, int v);
    void init_vector(double *vector, int NoOfElements);
    void init_integervector(int *vector, int dim);
    void init_matrix(double **matrix, int dim1, int dim2);
    void init_integermatrix(int **matrix, int dim1, int dim2);
    void init_3dimmatrix(double ***matrix, int dim1, int dim2, int dim3);
    void free_matrix(double **matrix, int dim1);
    void free_integermatrix(int **matrix, int dim1);
    void free_3dimmatrix(double ***matrix, int dim1, int dim2);
    void a_rungekutta(int p, double *avector, double **ka, double dt, double h, 
		  double **T);
    void rungekutta(int p, double *avector, double *gvector, double *bvector, 
		double **cmatrix, double dt, double h, double **T, double *t, 
		double **ka, double **kg, double **kb, double ***kc);
    double density(double t, int type, double *par);
    void ExportToMatlab_Phasetype(int dim, double h, double dt, double truncpoint,
			      double *pie, double **q, double *exit);
    void ExportToEMPHTmain_Phasetype(int dim, double h, double dt,double truncpoint,
				 double *pie, double **q, double *exit);
    void ExportToMatlab(int k, double *parameter, double truncpoint);
    void show_pi_T(int p, double *pi, double **T, double *t);
    double set_steplength(int p, double **T);
    void input_density();
    void input_sample(int NoOfInput[3]);
    void input_Csample(int NoOfInput[3]);
    void assign_vectors(int No);
    void assign_Cvectors(int NoOfInput[3]);
    int sort_observations(int size, double *vec1, double *vec2);
    int sort_interval(int NoOfPairs);
    int count_input();
    void Export_Sample();
    void Export_CensoredSample(int NoOfInput[3]);
    int AskForInput(int NoOfInput[3]);
    double rndom();
    void randomphase(int p, double *pi, double **T, double *exitvector, 
		 int *pilegal, int **Tlegal);
    void selectstructure(int p, double *pi, double **T, double *t, int *pilegal,
		     int **Tlegal);
    void EMstep(int p, double h, double *pi, double **T, double *t, double 
            *gvector, double *avector, double *bvector, double **cmatrix, 
            double *Bmean, double *Zmean, double **Nmean, double **kg, 
            double **ka, double **kb, double ***kc, double *ett, double 
            **g_left, double **a_left, double **b_left, double ***c_left);
    void compute_loglikelihood(double h, int p, double *pi, double **T, double *t, 
			   int stepindicator, double *avector, double **ka, 
			   double **a_left);
    void SavePhases(int p, double *pi, double **T);
    void EMiterate(int NoOfEMsteps, int p, double *pi, double **T, double *t,
	       double *gvector, double *avector, double *bvector, 
               double **cmatrix, double *Bmean, double *Zmean, double **Nmean, 
	       double **kg, double **ka, double **kb, double ***kc, 
               double *ett, double **g_left, double **a_left, 
               double **b_left, double ***c_left);
    int search_partner(void);

    
};

#endif /* PHASEFITTER_H */

