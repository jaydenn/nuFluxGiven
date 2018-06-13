#include <cmath>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "parameterStruct.h"
#include "nuRates.h"
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif

//natural log of Poisson dist: gives more accurate values for small probabilities (because of machine precision)
double logPoisson(double obs, double expect)
{
    if ( expect > 0. && obs > 0. )
        return -expect + obs * log( expect ) - gsl_sf_lngamma( obs+1 );
    else
        return -1E299;
}

//likelihood function for binned data
double logLikelihood(paramList *pList)
{

    //Calculate log-likelihood
    double loglike = 0;
    double Er_min, Er_max;
    double l,SM,BG;

    //loop over detectors
    for(int detj=0; detj < pList->ndet; detj++)
    {
        Er_min = pList->detectors[detj].ErL;
        //loop over recoil energy bins
        for(int i=0; i< pList->nBins; i++)
        {
            //set bin limits
            Er_max = Er_min + pList->detectors[detj].binW[i];
            SM  = intNuRate( Er_min, Er_max, pList, detj);
            BG  = pList->detectors[detj].BgNorm * intBgRate( &(pList->detectors[detj]), Er_min, Er_max);

            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(SM+BG) );
            loglike += l;
            //std::cout << "l: " << i << "/" << pList->nBins << " " << Er_min << " - " << Er_max << " expNu " << pList->detectors[detj].exposure*SM  << " expBG " << pList->detectors[detj].exposure*BG  << " obs: " << pList->detectors[detj].binnedData[i] << " " << l << std::endl;

            Er_min = Er_max; //update lower bin limit
        } 

    }
    if( std::isinf(loglike) )
        return -1e200;
    else
        return loglike;
}

void logLikelihoodGlobalFit(double *Cube, int &ndim, int &npars, double &lnew, long &pointer)    
{
    lnew=0;
    //get pointer in from MultiNest 
    paramList *pL = (paramList *) pointer;

    //scale pars for this point in the parameter space
    Cube[0] = pL->normPP = Cube[0]*3;
    Cube[1] = pL->normPEP= Cube[1]*15;
    Cube[2] = pL->normO  = Cube[2]*20;
    Cube[3] = pL->normN  = Cube[3]*20;
    Cube[4] = pL->normBE = Cube[4]*15;
//    Cube[5] = pL->normB  = Cube[5]*50;
//    Cube[6] = pL->normF  = Cube[6]*50;
    pL->normB =pL->normF = 1;

    //nuclear reaction chain priors
    if( pL->normPP + 2.36e-3*pL->normPEP < 8.49e-2*pL->normBE + 9.95e-5*pL->normB )
    {
        lnew=-1e299;
        return ;
    }
    if( 1.34*pL->normN < pL->normO )
    {
        lnew=-1e299;
        return ;
    }
    if( 37*pL->normO < pL->normF )
    {
        lnew=-1e299;
        return ;
    }


    lnew = - pow( pL->normPEP/pL->normPP - 1.006 ,2)/0.00034;

    //BG nuisance parameters
    int cubei=5;
    for(int detj=0; detj < pL->ndet; detj++)
    {
        if(pL->detectors[detj].BgUn > 1e-10)
        {
            pL->detectors[detj].BgNorm = 1;
            Cube[cubei] = pL->detectors[detj].krBgNorm = 1 + gsl_cdf_gaussian_Pinv( Cube[cubei], pL->detectors[detj].BgUn);
            Cube[cubei+1] = pL->detectors[detj].rnBgNorm = 1 + gsl_cdf_gaussian_Pinv( Cube[cubei+1], pL->detectors[detj].BgUn);
            Cube[cubei+2] = pL->detectors[detj].xeBgNorm = 1 + gsl_cdf_gaussian_Pinv( Cube[cubei+2], pL->detectors[detj].BgUn);
        }
    }

    //calculate total luminosity
    Cube[pL->nPar-1] = pL->normPP * 9.186e-1 + pL->normPEP*2.013e-3 + pL->normBE*7.388e-2 + pL->normN*1.2e-3 + pL->normO*5.641e-3 + pL->normF*1.53e-5 + pL->normB*4.339e-5;

    //calculate CNO fraction
    Cube[pL->nPar-2] = (pL->normN*1.2e-3 + pL->normO*5.641e-3 + pL->normF*1.53e-5)/(1.2e-3 + 5.641e-3 + 1.53e-5);

    //impose luminosity constraint?
    //if (pL->LC == 1)
    //   lnew += - pow( Cube[pL->nPar-1]-1,2) / (2.7e-7);
    
    lnew += logLikelihood(pL);

}


