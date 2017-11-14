#include <iostream>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>
#ifndef DETECTORSTRUCT_H
    #include "detectorStruct.h"
#endif    
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "nuRates.h"
#include "calcRates.h"

double detEff(double Er, int type)
{
    switch( type )
    {
        case 0: 
            return 1;
        case 1:
            return 0.331*(1+erf(0.248*(Er-9.22))); //SNS
        default:
            printf("invalid detector efficiency\n"); 
            return NAN;
    }
}

//return background in events/kg/day/keV
double detBackground(double Er, paramList *pList, int detj)
{
    
    switch( pList->detectors[detj].bg ) 
    {
        case 0: 
            return 1e-99;
        case 1: 
            return 0.01;
        case 2: 
            return 1e-5;
        default:
            printf("invalid detector background\n"); 
            return NAN; 
    }
}

double detRes(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 0;
        case 1: 
            return 0;
        case 2: 
            return 0;
        case 3: 
            return 0;
        case 4: 
            return 0;
        default:
            printf("invalid detector resolution\n"); 
            return NAN; 
    }
}

int newDetector(paramList *pList, char *name, double exp, int sourcej)
{
        if(pList->ndet==10)
        {
            std::cout << "already at max number of detectors (10)" << std::endl;
            return 1;
        }

        pList->detectors[pList->ndet].exposure = exp;
        pList->detectors[pList->ndet].sourcej = sourcej;
        
        //read in detector configuration
        FILE *detsINI;
        detsINI = fopen("detectors.ini","r");
        if(detsINI==NULL)
        {
            printf("unable to open detectors.ini\n");
            return 1;
        }

        char temp[200];
        char *ret;
        int err;

        ret = fgets(temp,200,detsINI);

        while(strcmp(temp,name)!=0)
        {
            err = fscanf(detsINI,"%s",temp);

            if(feof(detsINI))
            {
                printf("detector '%s' not found\n",name); 
                fclose(detsINI);
                return 1;
            }
        }

        sprintf( pList->detectors[pList->ndet].name, "%s", &(name[1]));

        while(temp[0]!='-')
        {
            err = fscanf(detsINI,"%s",temp);
 
            if(strcmp(temp,"Er")==0)  
                err=fscanf(detsINI,"%lf-%lf",&(pList->detectors[pList->ndet].ErL),&(pList->detectors[pList->ndet].ErU)); 
            if(strcmp(temp,"bg")==0)  
                err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].bg));
            if(strcmp(temp,"bgUn")==0)  
                err=fscanf(detsINI,"%lf",&(pList->detectors[pList->ndet].BgUn));
            if(strcmp(temp,"eff")==0) 
                err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].eff));            
            if(strcmp(temp,"res")==0) 
                err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].res));
        }

        ret = fgets(temp,200,detsINI);
        ret = fgets(temp,200,detsINI);
        ret = fgets(temp,200,detsINI);
        ret = fgets(temp,200,detsINI);

        while(!feof(detsINI) && temp[0]!='-')
        {    
            if(pList->detectors[pList->ndet].nIso==10)
            {
                std::cout << "already at max number of isotopes (10)" << std::endl;
                break;
            }
            sscanf(temp,"%d %d %lf %lf %lf %lf",&(pList->detectors[pList->ndet].isoZ[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoA[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoFrac[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoSZ[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoSN[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoJN[pList->detectors[pList->ndet].nIso])); 
            
            if ( pList->detectors[pList->ndet].isoJN[pList->detectors[pList->ndet].nIso] < 0.5 )
                pList->detectors[pList->ndet].isoJN[pList->detectors[pList->ndet].nIso] = 1e-99;
            
            pList->detectors[pList->ndet].nIso++;
            ret = fgets(temp,200,detsINI);           
        }
        ret = fgets(temp,200,detsINI);    
        ret = fgets(temp,200,detsINI);            

        std::string comma;
        int isoj=0;
        while( temp[0]!='_' && isoj < pList->detectors[pList->ndet].nIso)
        {
            int i = 0;
            std::istringstream ionizations(temp);

            while( std::getline(ionizations,comma,',') )
                pList->detectors[pList->ndet].ionization[isoj][i++] = atof(comma.c_str());

            isoj++;
               ret = fgets(temp,200,detsINI);    
        }   
        if(isoj > 1 && isoj != pList->detectors[pList->ndet].nIso)
        {
            std::cout << isoj << "Ionizations not listed for each isotope, aborting\n";
            return 1;
        }
        else if(isoj==1)
        {
            std::cout << "Ionizations will be duplicated for each isotope\n";
            
            for(int j=1; j<pList->detectors[pList->ndet].nIso; j++)
                 pList->detectors[pList->ndet].ionization[j] = pList->detectors[pList->ndet].ionization[0];
        }
        
        //finished reading in det data         
        fclose(detsINI);
        
        //only need to calculate background and SM signal once, store in a table for interpolation, stored as events/kg/day/keV
        //get values of bg at relevant energies
        std::cout << "Initializing rates\n";
        nuRatesInit( pList, pList->ndet);
        rateInit( pList, pList->ndet, &detBackground, pList->detectors[pList->ndet].background);

        //must initialize each flux in turn
        //for(int fluxj=0; fluxj< pList->sources[pList->detectors[pList->ndet].sourcej].numFlux; fluxj++)
        
        std::cout << "\ndone." << std::endl; 
        
        pList->ndet++;
        
        return 0;
        
}

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector *det, double Er_min, double Er_max)                            
{   
    if( Er_min > det->ErL && Er_max < det->ErU )
        return gsl_spline_eval_integ(det->background, Er_min, Er_max, det->accelBg);
    else if( Er_min < det->ErL && Er_max < det->ErU )
        return gsl_spline_eval_integ(det->background, det->ErL, Er_max, det->accelBg);
    else
        return gsl_spline_eval_integ(det->background, det->ErL, det->ErU, det->accelBg);
}

double diffBgRate(detector det, double Er)
{
    return gsl_spline_eval(det.background, Er, det.accelBg);
}
