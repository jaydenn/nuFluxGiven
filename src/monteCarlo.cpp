#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "nuRates.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif

int SEED=0;

void generateBinnedData(paramList *pList, int detj, int b, int simSeed)
{

    double Er_min, Er_max;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);
    pList->detectors[detj].nEvents = 0;

    try
    {
        pList->detectors[detj].binnedData = new double[pList->nBins];
        pList->detectors[detj].binW = new double[pList->nBins];
    }
    catch (std::bad_alloc& ba)
    {
      std::cout << "bad_alloc caught: " << ba.what() << std::endl << "you requested: " << pList->nBins << " doubles" <<std::endl;
      return ;
    }

    Er_min = pList->detectors[detj].ErL;
    double logBinW = ( log10( pList->detectors[detj].ErU ) - log10 (pList->detectors[detj].ErL ) ) / ( (double) pList->nBins);

    double ranBG,ranFlux,BG,SM;
    if( pList->asimov == 0 )
    {
        ranFlux= 1+gsl_ran_gaussian(r, pList->sources[pList->detectors[pList->detj].sourcej].nuFluxUn[0]);
        ranBG  = 1+gsl_ran_gaussian(r, pList->detectors[detj].BgUn);
    }

    for(int i=0; i<pList->nBins; i++)
    {

        if(pList->logBins==1)
            pList->detectors[detj].binW[i] = pow(10, log10(Er_min) + logBinW) - Er_min;
        else
            pList->detectors[detj].binW[i] = ( pList->detectors[detj].ErU - pList->detectors[detj].ErL ) / ( (double) pList->nBins);

        Er_max = Er_min + pList->detectors[detj].binW[i];

        SM  = intNuRate( Er_min, Er_max, pList, detj);
        BG  = b * intBgRate(&(pList->detectors[detj]), Er_min, Er_max) ;

        if( pList->asimov == 1)
            pList->detectors[detj].binnedData[i] = pList->detectors[detj].exposure * ( SM + BG );
        else
            pList->detectors[detj].binnedData[i] = gsl_ran_poisson(r, pList->detectors[detj].exposure *( ranFlux*SM + ranBG*BG ));
        //std::cout << " MC" << i << " " << Er_min << " - " << Er_max << ": bg " << BG <<  " sm " << SM << " obs " << pList->detectors[detj].binnedData[i] << std::endl;
        pList->detectors[detj].nEvents += pList->detectors[detj].binnedData[i];

        Er_min = Er_max; //update lower bin limit

    }

    return ;

}
