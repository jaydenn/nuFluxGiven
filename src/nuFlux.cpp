#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <string>
#include <fstream>
#include <sstream>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif  
#ifndef SOURCESTRUCT_H
    #include "sourceStruct.h"
#endif  

const double HBARC = 1.975e-14; //GeV*cm
const double GEVMETER = HBARC/100.0; //GeVcm
double FIRSTEVALCNS = 0;

gsl_integration_workspace * W;

int sourceInit(paramList *pL, char *sourceName, double sourceDistance)
{
    //check if source already exists
    int sourcej=0;
    while(sourcej < pL->nSource)
    {
        if(strcmp(sourceName,pL->sources[sourcej].name)==0)
            return sourcej;
        sourcej++;
    }
    //if not create new one
    pL->nSource++;
    strcpy(pL->sources[sourcej].name,sourceName);
    //read in source file
    std::ifstream sourceFile("sources.ini");
    
    std::string line;
    std::getline(sourceFile, line);
    
    while( line.compare(sourceName) != 0)
    {
        if(sourceFile.eof())
        {    
            std::cout << "source " << sourceName << " not found\n";
            return -1; 
        } 
        else
            std::getline(sourceFile, line);
    }
    
    //read in component fluxes
    std::string fluxFile;
    int fluxj=0;
    double fluxE,fluxN,lineEnergy;
    std::string plusMinus = "+/-";
    std::string flavor;
    std::getline(sourceFile, line);
    while( line[0] != '-' && fluxj < 10 && !sourceFile.eof())
    {
        
        std::istringstream lineStream(line);
        //is this a lineFlux (in energy)?        
        if(line.compare(0,4,"line")==0)
        {
            
            if(!( lineStream >> fluxFile >> lineEnergy >> flavor >> pL->sources[sourcej].nuFlux[fluxj] >> plusMinus >> pL->sources[sourcej].nuFluxUn[fluxj] ))
            {
                std::cout << "error parsing source data (line)\n";
                return -1;
            }
            pL->sources[sourcej].nuFluxUn[fluxj] /= pL->sources[sourcej].nuFlux[fluxj]; //want fractional uncertainty
            pL->sources[sourcej].nuFlux[fluxj] /= pow(sourceDistance,2);
            pL->sources[sourcej].nuFlux[fluxj] *= pow(HBARC,2);               //convert units
            pL->sources[sourcej].lineE[fluxj] = lineEnergy/1000;
            pL->sources[sourcej].isLine[fluxj] = 1;
            
        }
        else
        {

            if(!( lineStream >> fluxFile >> flavor >> pL->sources[sourcej].nuFlux[fluxj] >> plusMinus >> pL->sources[sourcej].nuFluxUn[fluxj] ))
            {
                std::cout << " error parsing source data\n";
                return -1;
            }
            
            pL->sources[sourcej].nuFluxUn[fluxj] /= pL->sources[sourcej].nuFlux[fluxj];  //want fractional uncertainty
            pL->sources[sourcej].nuFlux[fluxj] /= pow(sourceDistance,2);   
            pL->sources[sourcej].nuFlux[fluxj] *= pow(HBARC,2);           //convert units from /cm^2
            
            if( !( lineStream >> pL->sources[sourcej].survProb[fluxj] >> plusMinus >> pL->sources[sourcej].survProbUn[fluxj] ) )
            {
                pL->sources[sourcej].survProb[fluxj]=1;
                pL->sources[sourcej].isSolar[fluxj]=0;
            }
            else
            {
                pL->sources[sourcej].isSolar[fluxj]=1;
                pL->sources[sourcej].nuFluxUn[fluxj] = sqrt( pow(pL->sources[sourcej].nuFluxUn[fluxj],2) + pow(pL->sources[sourcej].survProbUn[fluxj]/pL->sources[sourcej].survProb[fluxj],2) );
            }   
            
            std::cout<< "reading " << fluxFile << std::endl;
            std::ifstream flux(fluxFile.c_str());
            double *fluxP_E = new double [500];
            double *fluxP_N = new double [500];
            
            int i=0;
            int maxPoints=500;
            while( flux >> fluxE >> fluxN )
            {
                fluxP_E[i  ] = fluxE/1000;   //convert to GeV
                fluxP_N[i++] = fluxN*1000;
                pL->sources[sourcej].flux_points[fluxj]++;     

                if(i==maxPoints)
                {
                    maxPoints+=100;
                    fluxP_E = (double *) realloc(fluxP_E, maxPoints * sizeof(double));
                    fluxP_N = (double *) realloc(fluxP_N, maxPoints * sizeof(double));
                }       
            }
            flux.close(); flux.clear();
            //initialize gsl interpolator for flux
            
            pL->sources[sourcej].nuFluxInterp[fluxj] = gsl_spline_alloc(gsl_interp_linear, pL->sources[sourcej].flux_points[fluxj]);
            pL->sources[sourcej].nuFluxAccel[fluxj]  = gsl_interp_accel_alloc();
            
            gsl_spline_init(pL->sources[sourcej].nuFluxInterp[fluxj], fluxP_E, fluxP_N, pL->sources[sourcej].flux_points[fluxj]);
            pL->sources[sourcej].EnuMin[fluxj] = fluxP_E[0];
            pL->sources[sourcej].EnuMax[fluxj] = fluxP_E[pL->sources[sourcej].flux_points[fluxj]-1];
            
            double norm = gsl_spline_eval_integ (pL->sources[sourcej].nuFluxInterp[fluxj], pL->sources[sourcej].EnuMin[fluxj], pL->sources[sourcej].EnuMax[fluxj], pL->sources[sourcej].nuFluxAccel[fluxj]);

            if ( fabs(norm-1.00) > .01)
            {    
                std::cout << "ERROR: flux data in " << fluxFile << " is not properly normalized N = " << norm << std::endl;
                return -1;
            }
        }
        
        if( flavor == "e" )
            pL->sources[sourcej].nuFluxFlav[fluxj] = 1;
        else if( flavor == "ebar")
            pL->sources[sourcej].nuFluxFlav[fluxj] = -1;
        else if( flavor == "mu" )
            pL->sources[sourcej].nuFluxFlav[fluxj] = 2;
        else if( flavor == "mubar" )
            pL->sources[sourcej].nuFluxFlav[fluxj] = -2;
        else if( flavor == "tau" )
            pL->sources[sourcej].nuFluxFlav[fluxj] = 3;
        else if( flavor == "taubar" )
            pL->sources[sourcej].nuFluxFlav[fluxj] = -3;
        else
        {
            std::cout << "flavor of neutrino flux not recognized\n";
            return -1;
        }
                    
        fluxj++;
        std::getline(sourceFile, line);

    }
    if(fluxj==10)
        std::cout << "WARNING: max number of flux elements reached (10), ignoring any further components\n";
    
    pL->sources[sourcej].numFlux = fluxj;
    sourceFile.close();

    return sourcej;
    
}

//returns diffNuFlux in GeV per sec, at the point Enu(GeV) ( /cm^2/s/GeV * hc^2)
double nuFlux(double EnuGeV, paramList *pL, int sourcej, int fluxj)
{
    if(EnuGeV < pL->sources[sourcej].EnuMax[fluxj] && EnuGeV > pL->sources[sourcej].EnuMin[fluxj])
        return pL->sources[sourcej].nuFlux[fluxj] * gsl_spline_eval(pL->sources[sourcej].nuFluxInterp[fluxj], EnuGeV, pL->sources[sourcej].nuFluxAccel[fluxj]);
    else
        return 1e-299;
}


//units GeV/s
double EnuIntegrand0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj);
}

//units 1/s
double EnuIntegrand1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj) / EnuGeV ;
}

//units 1/GeV/s
double EnuIntegrand2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj) / ( EnuGeV*EnuGeV ) ;
}

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow, int sourcej, int fluxj)
{
	int limit = 3000;
	double integral,absErr,tol;
	
	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (3000);
		FIRSTEVALCNS=1;
	}
	
	if(EnuPow==0)
	{
		pList->F.function = &EnuIntegrand0;
		tol=1e-21;
	}
	else if(EnuPow==-1)
	{
		pList->F.function = &EnuIntegrand1;
		tol=1e-20;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrand2;
		tol=1e-19;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
	
	double EnuMinGeV = 0.5 * (ErGeV + sqrt( pow(ErGeV,2) + 2*ErGeV*Mt ) );
    
    pList->sourcej = sourcej;
    pList->fluxj = fluxj;

    if(pList->sources[sourcej].isLine[fluxj]==1)
    {
        if(EnuMinGeV < pList->sources[sourcej].lineE[fluxj] )
            integral = pList->sources[sourcej].nuFlux[fluxj] * pow( pList->sources[sourcej].lineE[fluxj], EnuPow);
        else
            integral = 1e-199;
    }
    else
        gsl_integration_qag(&(pList->F), EnuMinGeV, pList->sources[sourcej].EnuMax[fluxj], tol, 1e-3, limit, 1, W, &integral, &absErr); 

    if (integral < 0)        
        return 1e-299;
    else    
    	return pList->sources[sourcej].nuFluxNorm[fluxj]*integral;
		
}


//units GeV/s
double EnuIntegrandOsc0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1.0-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->sources[pList->sourcej].distance/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj);
}

//units 1/s
double EnuIntegrandOsc1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1.0-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->sources[pList->sourcej].distance/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj) / EnuGeV ;
}

//units 1/GeV/s
double EnuIntegrandOsc2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->sources[pList->sourcej].distance/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->sourcej, pList->fluxj) / ( EnuGeV*EnuGeV ) ;
}

//se
double fluxIntegralOsc(double ErGeV,  paramList *pList, double Mt, int EnuPow, int sourcej, int fluxj)
{
	int limit = 3000;
	double integral,absErr,tol;
	
	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (3000);
		FIRSTEVALCNS=1;
	}
	
	if(EnuPow==0)
	{
		pList->F.function = &EnuIntegrandOsc0;
		tol=1e-21;
	}
	else if(EnuPow==-1)
	{
		pList->F.function = &EnuIntegrandOsc1;
		tol=1e-19;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrandOsc2;
		tol=1e-18;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
	
	double EnuMinGeV = 0.5 * (ErGeV + sqrt( pow(ErGeV,2) + 2*ErGeV*Mt ) );
    
    pList->sourcej = sourcej;
    pList->fluxj = fluxj;
        
    if(pList->sources[sourcej].isLine[fluxj]==1)
    {
        if(EnuMinGeV < pList->sources[sourcej].lineE[fluxj] )
            integral = pList->sources[sourcej].nuFlux[fluxj] * pow( pList->sources[sourcej].lineE[fluxj], EnuPow);
        else
            integral = 1e-199;
    }
    else
        gsl_integration_qag(&(pList->F), EnuMinGeV, pList->sources[sourcej].EnuMax[fluxj], tol, 1e-3, limit, 1, W, &integral, &absErr); 

    if (integral < 0)        
        return 1e-299;
    else    
    	return pList->sources[sourcej].nuFluxNorm[fluxj]*integral;
		
}

