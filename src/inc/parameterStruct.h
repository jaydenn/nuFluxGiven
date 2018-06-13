//definition of parameter structs
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#ifndef IOSTREAM
	#include <iostream>
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#ifndef SOURCESTRUCT_H
    #include "sourceStruct.h"
#endif  
#define PARAMETERSTRUCT_H

struct paramList {

    //root directory for file output
    char root[50];
    int nPar;
    //MultiNest sampling parameters
    double sampling[8];
    //Struct for neutrino sources
    int nSource, sourcej, fluxj;
    sourceStruct sources[10];
    double Er;
    int LC;
    int nucPriors;

	//setup for rate integration
	gsl_function F;
	double (*rateFunc)(double, double, paramList *, int);

    double normPP;
    double normPEP;
    double normO;
    double normN;
    double normBE;
    double normB;
    double normF;

    //BSM model parameters
    int nucScat;
    int elecScat;

    double maxL;
    int asimov;
    int logBins;
    int maxBins;
    int nBins;
    
	int ndet, detj;
	detector detectors[10];

	void printPars()
	{
   		std::cout << "Conudl configuration:" << std::endl;		
    	std::cout << "  root: " << root << std::endl;
    	std::cout << "  asimov: " << asimov << std::endl;
    	std::cout << "  logBins: " << logBins << std::endl;
    	//std::cout << "  flux: " << nuFlux << " +/- " << nuFluxUn*100 << "%" << std::endl;
		std::cout << ndet << " Detector(s):" << std::endl;
		for(int i=0;i<ndet;i++)
			detectors[i].printDetSpecs();
	    std::cout << ndet << " Source(s):" << std::endl;
	   // for(int i=0;i<nSource;i++)
	//	sources[i].printSourceSpecs();
    }

	paramList()
	{
    	normPP=normPEP=normO=normN=normBE=normB=normF=1;    	    
	    elecScat=nucScat=0;
		ndet=detj=fluxj=nSource=sourcej=LC=0;
        asimov=logBins=1;
        maxBins = 100;
	}
	
};
	
