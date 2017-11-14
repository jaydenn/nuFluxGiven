#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <gsl/gsl_interp.h>
#ifndef GSL_SPLINE_H
    #include <gsl/gsl_spline.h>
#endif
#include "detectorStruct.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
int nuRatesInit( paramList *pL, int detj)
{

    //read in rate files
    char const *files[] = {"data/ppn_xe.dat","data/pep_xe.dat","data/N13_xe.dat","data/O15_xe.dat","data/Be7_xe.dat","data/B8_xe.dat","data/F17_xe.dat"};
    gsl_spline *splines[] = {pL->detectors[detj].NR_PP, pL->detectors[detj].NR_PEP, pL->detectors[detj].NR_N, pL->detectors[detj].NR_O, pL->detectors[detj].NR_BE,pL->detectors[detj].NR_B,pL->detectors[detj].NR_F};
    std::ifstream RFF;

    double val;
 
    double Er[10000];
    double Rate[10000];

    for (int j=0; j<7; j++)
    {
        RFF.open(files[j]);
        std::cout << "reading " << files[j] << std::endl;
        
        if(!RFF)
        {
            std::cout << "file not found for detector " << pL->detectors[detj].name << "\n";
            return 1;
        }   
        for (int i=0; i<10000; i++)
        {
            RFF >> Er[i] >> Rate[i];
            Rate[i]*=0.03337;
            Er[i]*=1000;
        }
        
        gsl_spline_init(splines[j],Er,Rate,10000);
        RFF.close(); RFF.clear();
    }            
    
    return 0;
    
}

double diffNuRate(double Er, paramList *pL, int detj) 
{
    double rate = 1e-99;
    rate += pL->normPP  * gsl_spline_eval(pL->detectors[detj].NR_PP, Er, pL->detectors[detj].accelPP);
    rate += pL->normPEP * gsl_spline_eval(pL->detectors[detj].NR_PEP,Er, pL->detectors[detj].accelPEP);
    rate += pL->normO   * gsl_spline_eval(pL->detectors[detj].NR_O,  Er, pL->detectors[detj].accelO);
    rate += pL->normN   * gsl_spline_eval(pL->detectors[detj].NR_N,  Er, pL->detectors[detj].accelN);
    rate += pL->normBE  * gsl_spline_eval(pL->detectors[detj].NR_BE, Er, pL->detectors[detj].accelBE);
    rate += pL->normB  * gsl_spline_eval(pL->detectors[detj].NR_B, Er, pL->detectors[detj].accelB);
    rate += pL->normF  * gsl_spline_eval(pL->detectors[detj].NR_F, Er, pL->detectors[detj].accelF);

    return rate;
}

//placeholder for future interpolating function if required
double intNuRate(double Er_min, double Er_max, paramList *pL, int detj) 
{
    double rate = 1e-99;
    rate += pL->normPP  * gsl_spline_eval_integ(pL->detectors[detj].NR_PP, Er_min, Er_max, pL->detectors[detj].accelPP);
    rate += pL->normPEP * gsl_spline_eval_integ(pL->detectors[detj].NR_PEP,Er_min, Er_max, pL->detectors[detj].accelPEP);
    rate += pL->normO   * gsl_spline_eval_integ(pL->detectors[detj].NR_O,  Er_min, Er_max, pL->detectors[detj].accelO);
    rate += pL->normN   * gsl_spline_eval_integ(pL->detectors[detj].NR_N,  Er_min, Er_max, pL->detectors[detj].accelN);
    rate += pL->normBE  * gsl_spline_eval_integ(pL->detectors[detj].NR_BE, Er_min, Er_max, pL->detectors[detj].accelBE);
    rate += pL->normB  * gsl_spline_eval_integ(pL->detectors[detj].NR_B, Er_min, Er_max, pL->detectors[detj].accelB);
    rate += pL->normF  * gsl_spline_eval_integ(pL->detectors[detj].NR_F, Er_min, Er_max, pL->detectors[detj].accelF);
    
    return rate;
}

