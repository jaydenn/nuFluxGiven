#include <iostream>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
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

//return background in events/t/year/keV
double detBackground(double Er, paramList *pList, int detj)
{
    
    switch( pList->detectors[detj].bg ) 
    {
        case 0: 
            return 1e-99;
        case 1:
            return 55;
        case 2: 
            return 1e-5;
        case 3:
            return gsl_spline_eval(pList->detectors[detj].background, Er, pList->detectors[detj].accelBg);
        default:
        {
            printf("invalid detector background\n"); 
            return NAN; 
        }
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

        if(pList->detectors[pList->ndet].bg == 3)
        {
            std::ifstream RFF;
           
            RFF.open("data/2nuBB_smear10percent.dat");
            double Er[2400];
            double Rate[2400];
            for (int i=0; i<2400; i++)
            {    
                RFF >> Er[i] >> Rate[i];
                Rate[i]*=0.001;
            }
            pList->detectors[pList->ndet].xeBackground = gsl_spline_alloc(gsl_interp_linear,2400);
            gsl_spline_init( pList->detectors[pList->ndet].xeBackground,Er,Rate,2400);
            RFF.close();
            
            std::ifstream KR;
            KR.open("data/Kr_smear10percent.dat");
            for (int i=0; i<620; i++)
                KR >> Er[i] >> Rate[i];
            pList->detectors[pList->ndet].krBackground = gsl_spline_alloc(gsl_interp_linear,620);
            gsl_spline_init( pList->detectors[pList->ndet].krBackground,Er,Rate,620);
            KR.close();
            
            std::ifstream RN;
            RN.open("data/Rn_smear10percent.dat");
            for (int i=0; i<947; i++)
                RN >> Er[i] >> Rate[i];
            pList->detectors[pList->ndet].rnBackground = gsl_spline_alloc(gsl_interp_linear,947);
            gsl_spline_init( pList->detectors[pList->ndet].rnBackground,Er,Rate,947);
            RN.close();
            
        }
        else
        {
            pList->detectors[pList->ndet].background = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
            rateInit( pList, pList->ndet, &detBackground, pList->detectors[pList->ndet].background);
        }
        //must initialize each flux in turn
        //for(int fluxj=0; fluxj< pList->sources[pList->detectors[pList->ndet].sourcej].numFlux; fluxj++)

        std::cout << "\ndone." << std::endl; 

        pList->ndet++;

        return 0;
        
}

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector *det, double Er_min, double Er_max)                            
{   
    if(det->bg == 3)
    {
        double totBG = 0;
        if( Er_min >= det->ErL && Er_max <= 620 )
            totBG+= det->krBgNorm*gsl_spline_eval_integ(det->krBackground, Er_min, Er_max, det->accelkrBg);
        else if( Er_min < det->ErL && Er_max < det->ErU )
            totBG+=  det->krBgNorm*gsl_spline_eval_integ(det->krBackground, det->ErL, Er_max, det->accelkrBg);
        else if( Er_min < det->ErL && Er_max > 620 )
            totBG+=  det->krBgNorm*gsl_spline_eval_integ(det->krBackground, det->ErL, det->ErU, det->accelkrBg);
        if( Er_min >= det->ErL && Er_max <= 950 )
            totBG+=  det->rnBgNorm*gsl_spline_eval_integ(det->rnBackground, Er_min, Er_max, det->accelrnBg);
        else if( Er_min < det->ErL && Er_max < 950 )
            totBG+=  det->rnBgNorm*gsl_spline_eval_integ(det->rnBackground, det->ErL, Er_max, det->accelrnBg);
        else if( Er_min < det->ErL && Er_max > 950 )
            totBG+=  det->rnBgNorm*gsl_spline_eval_integ(det->rnBackground, det->ErL, 950, det->accelrnBg);
        if( Er_min >= det->ErL && Er_max <= det->ErU )
            totBG+=  det->xeBgNorm*gsl_spline_eval_integ(det->xeBackground, Er_min, Er_max, det->accelxeBg);
        else if( Er_min < det->ErL && Er_max < det->ErU )
            totBG+=  det->xeBgNorm*gsl_spline_eval_integ(det->xeBackground, det->ErL, Er_max, det->accelxeBg);
        else if( Er_min < det->ErL && Er_max > det->ErU )
            totBG+=  det->xeBgNorm*gsl_spline_eval_integ(det->xeBackground, det->ErL, det->ErU, det->accelxeBg);

        return totBG;
    }
    else
    {
        if( Er_min >= det->ErL && Er_max <= det->ErU )
            return det->BgNorm*gsl_spline_eval_integ(det->background, Er_min, Er_max, det->accelBg);
        else if( Er_min < det->ErL && Er_max < det->ErU )
            return det->BgNorm*gsl_spline_eval_integ(det->background, det->ErL, Er_max, det->accelBg);
        else if( Er_min < det->ErL && Er_max > det->ErU )
            return det->BgNorm*gsl_spline_eval_integ(det->background, det->ErL, det->ErU, det->accelBg);
        else 
            return 0;
    }
}

double diffBgRate(detector det, double Er)
{
    if( det.bg=3 )
    {
        double totBG=0;
        if( Er > 5 && Er < 620 )
            totBG+=gsl_spline_eval(det.krBackground, Er, det.accelkrBg);
        if( Er > 5 && Er < 950 )
            totBG+=gsl_spline_eval(det.rnBackground, Er, det.accelrnBg);
        if( Er > 5 && Er < 2400 )
            totBG+=gsl_spline_eval(det.xeBackground, Er, det.accelxeBg);
        return totBG;
    }
    else
    {
        return gsl_spline_eval(det.background, Er, det.accelBg);
    }
}
