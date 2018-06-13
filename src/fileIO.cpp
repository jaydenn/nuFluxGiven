#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#ifndef SOURCESTRUCT_H
	#include "sourceStruct.h"
#endif	

//gets sampling parameters from file
int readConfigFile(paramList *pL, char *filename) 
{

    FILE* input;
    input = fopen(filename,"r");
    if(input==NULL) 
    {
        std::cout << "unable to open parameter file: " << filename << std::endl;
        return -1;	
    }
  
    int mode;
    char *ret;
    char temp[400];
    char root[50];
    ret = fgets(temp,200,input);

    //Mode switch
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&mode);
    
    if ( mode == 2 )
    {
        //Get Multinest sampling parameters
        FILE* input;
        input = fopen("multinest.ini","r");
        if(input==NULL) 
        {
            std::cout << "unable to open multinest.ini\n";
            return -1;	
        }
        ret = fgets(temp,200,input);
        for (int i=0;i<8;i++)
        {
            ret = fgets(temp,200,input);
            sscanf(temp,"%lf",&(pL->sampling[i]));
        }
    }
    
    // root for output files
    ret = fgets(temp,200,input);
    sscanf(temp,"%s %*s",root);
    sprintf(pL->root, "%s", root);		 
  
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->nBins));

    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->logBins));
    
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->nucPriors));

    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->LC));
    
    //Detector setup
    double exp;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);

    char name[30];
    char sourceName[30]="";
    int sourcej;
    double sourceDistance=0;
    
    while(temp[0]=='#')
    {
        sscanf(temp,"%s %lf %s %lf", name, &exp, sourceName, &(sourceDistance));
        //sourcej=sourceInit(pL, sourceName, sourceDistance);
        //if(sourcej<0) { std::cout << "Could not create source, exiting" << std::endl; return -1; }
        
	    if(newDetector(pL, name, exp, sourcej)) { std::cout << "Could not create detector, exiting" << std::endl; return -1; }
        std::cout << "Using detector " << name << " with source " << sourceName << std::endl;
        ret = fgets(temp,200,input);
    }    
    
    //asimov or random sim?
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->asimov));
    
    fclose(input); 

    return mode;
}

