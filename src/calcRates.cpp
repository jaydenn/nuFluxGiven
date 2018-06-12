#include <iostream>
#include <fstream>
#include <iomanip>
#include "nuRates.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int calcRates(paramList *pList)
{
    double ErkeV;
    char filename[90];    
    
    //format output streams
    std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    
    for(int detj=0; detj < pList->ndet; detj++)
    {
        std::cout << "------------------------\n";
        std::cout << detj+1 << ". " << pList->detectors[detj].name << std::endl;
        std::cout << "------------------------\n";
        std::cout << "  total rates: \n"; 
        std::cout << "     SM  = " << intNuRate(  pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj)         << " events/kg/day" << std::endl;
        std::cout << "     BG  = " << intBgRate(  &(pList->detectors[detj]), pList->detectors[detj].ErL, pList->detectors[0].ErU) << " events/kg/day" << std::endl;
        
        std::cout << "  differential rates: \n"; 
        std::cout << "    Er (keV)        SM dN/dE        BG dN/dE    (events/kg/day/keV)" << std::endl;

        for (int i=0; i<501; i+=1)
        {
            if(pList->logBins == 1)
                ErkeV = pow(10, log10(pList->detectors[detj].ErL) + (double)i*(log10(pList->detectors[detj].ErU)-log10(pList->detectors[detj].ErL))/500);
            else
                ErkeV = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/500;

            std::cout << "    " << ErkeV << "      " << diffNuRate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV) << std::endl;
        }
        
    }
    
}    

void rateInit( paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline)
{
    double ErkeV[INTERP_POINTS];
    double rate[INTERP_POINTS];

    double linStep,logStep; 
    
    logStep = pow(pList->detectors[detj].ErU/pList->detectors[detj].ErL/0.98,1/(INTERP_POINTS-10.0));
    linStep = (pList->detectors[detj].ErU-pList->detectors[detj].ErL*0.98)/(INTERP_POINTS-10.0);

    ErkeV[0] = 0.99*pList->detectors[detj].ErL; 
    rate[0] = rateFunc( (double)ErkeV[0], pList, detj);
    
    for( int i=1; i < INTERP_POINTS; i++ )
    {
        //always over and undershoot range so that interpolation is well behaved
        if(pList->logBins == 0 )//&& ErkeV[i-1] > 5)
            ErkeV[i] = ErkeV[i-1] + linStep;
        else
            ErkeV[i] = ErkeV[i-1] * logStep;
            
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj);	
        //std::cout << i << " " << ErkeV[i] << " " << logStep << std::endl;
    }
    if(ErkeV[INTERP_POINTS-1] < pList->detectors[detj].ErU)
    {
        std::cout << ErkeV[INTERP_POINTS-1] << " < " << pList->detectors[detj].ErU << " interpolation failed to cover range adequately, please check\n";
        ErkeV[INTERP_POINTS-1] = 1.01*pList->detectors[detj].ErU;
        rate[INTERP_POINTS-1] = rateFunc( (double)ErkeV[INTERP_POINTS-1], pList, detj);
    }
    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,INTERP_POINTS);
}


