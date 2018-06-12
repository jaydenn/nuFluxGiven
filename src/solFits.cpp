#include <iostream>
#include <cstring>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "fileIO.h"
#include "calcRates.h"
#include "globalFit.h"


int main(int argc, char *argv[])
{
        
    //argument validation and file input
    char filename[100];
    int err=0;
    int mode=0;
    paramList pList;
    
    if(argc==1)
    {
        std::cout << "reading from default config.dat file" << std::endl; 
        sprintf(filename,"config.dat"); 
        mode = readConfigFile(&pList,filename);
    }
    else if(argc>1)
    {
        for(int i=1; i<argc; i++ )
        {
            
            if ( std::strcmp( argv[i], "-c") == 0)
            {
               sprintf(filename,"%s",argv[i+1]);
               mode = readConfigFile(&pList,filename);
               i++;
            }
            else
            {
                std::cerr << "Invalid input" << std::endl << std::endl;
                std::cerr << "Usage: ./solFits " << std::endl;
                std::cerr << "       default (no flags):" << std::endl; 
                std::cerr << "               (solFits runs with the default \"config.dat\" parameter file)" << std::endl;
                std::cerr << "       optional flags: " << std::endl; 
                std::cerr << "               -c file.dat (solFits    runs with the specified parameter file)" << std::endl << std::endl;
                break;
            }
        }
    } 
    
	//pList.printPars();
    
	if ( mode < 1 ) 
    {
        std::cerr << "Problem with configuration file, aborting" << std::endl;
        return 0;
    }
	
	//print-out rate mode
	if(mode == 1)
    {
        calcRates(&pList);
        return 0;
    }   

    //NSI multinest fit
    if (mode == 2 )
    {
        globalFit(&pList);
        return 0;
    }
    
    std::cerr << "choose a valid mode" << std::endl;
    return 0;
}
