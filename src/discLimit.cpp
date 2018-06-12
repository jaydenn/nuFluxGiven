#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include <string.h>
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif
#include "likelihood.h"
#include "monteCarlo.h"
#include "nuRates.h"

double my_LS(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->normO  = pL->normN  = gsl_vector_get(v, 0);
    
    pL->detectors[pL->detj].BgNorm = fabs(gsl_vector_get(v, 1));
    l += 0.5/pow(pL->detectors[0].BgUn,2) * pow(pL->detectors[0].BgNorm - 1,2);
    
    pL->normPP = fabs(gsl_vector_get(v, 2));
    pL->normPEP= fabs(gsl_vector_get(v, 3));
    pL->normBE = fabs(gsl_vector_get(v, 4));
    
    l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normPP - 1,2) \
      +  0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normPEP - 1,2) \
      +  0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normBE - 1,2);
    
    l -= logLikelihood(pL);
    //std::cout << " like "  << l << " " << pL->signalNorm << " " << pL->source.nuFluxNorm[0]  << "\n";
    return l;

}

double my_L0(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->detectors[pL->detj].BgNorm  = gsl_vector_get(v, 0);
    l += 0.5/pow(pL->detectors[0].BgUn,2) * pow(pL->detectors[0].BgNorm - 1,2);
    
    pL->normPP = fabs(gsl_vector_get(v, 1));
    pL->normPEP= fabs(gsl_vector_get(v, 2));
    pL->normBE = fabs(gsl_vector_get(v, 3));
    
    l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normPP - 1,2) \
      +  0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normPEP - 1,2) \
      +  0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->normBE - 1,2);
    
    l -= logLikelihood(pL);
    //std::cout << l << std::endl;
    return l;

}

double findMaxLS(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = my_LS;
    my_func.params = (void *)pL;

    //start point and step size
    gsl_vector *x,*dx;
    my_func.n = 5;

    x = gsl_vector_alloc ( my_func.n );
    dx = gsl_vector_alloc ( my_func.n );
    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/10);
    
    for(int i=1; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .01);
    }
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        
        //std::cout << "       " <<iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " << gsl_vector_get (s->x, 2) << " L " << s->fval << std::endl; 
    }
    while (iter < 1000 && gsl_multimin_fminimizer_size(s)/s->fval > .00001);
    
   if(iter==1000)
        std::cout << "LS non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval << " > .0005  " << std::endl;
    
    double LS =  s->fval;
    pL->signalNorm = gsl_vector_get(s->x, 1);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return -LS;
}


double findMaxL0(paramList *pL)
{

    size_t iter = 0;
    int status;
    
    double mu = pL->signalNorm; //save current mu for later
    pL->signalNorm = 0;
    
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    my_func.f = my_L0;
    my_func.params = (void *)pL;

    //start point
    gsl_vector *x,*dx;
    my_func.n = 1 + pL->source.numFlux;
    
    x = gsl_vector_alloc (my_func.n);
    dx = gsl_vector_alloc (my_func.n);
    
    for(int i=0; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .0001);
    }
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
       // std::cout << "       " << iter << " " <<  gsl_vector_get (s->x, 0) << " " << gsl_vector_get (s->x, 1) << " " << s->fval << std::endl; 
    }
    while (iter < 600 && gsl_multimin_fminimizer_size(s)/s->fval > 1e-9 && !status);
    if(iter==600)
        std::cout << "L0 non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval  << " > 1e-9 " <<  std::endl;

    double L0 = s->fval;
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    pL->signalNorm = mu;
    
    return -L0;
}

//statistic for discovery
double q0(paramList *pL)
{	 

    double maxL0 = findMaxL0( pL );
    double maxL;
    //if using asimov just set parameters to MLE
    if(pL->asimov==1)
    {
        pL->normO = pL->normN = 1;
        pL->detectors[pL->detj].BgNorm = 1;
        pL->normPP = 1;
        pL->normPEP= 1;
        pL->normBE = 1;
        maxL = logLikelihood(pL);
    }
    else
        maxL = findMaxLS( pL );
    
    double q = - 2 * ( maxL0 - maxL ); 
    //std::cout << "returning q= " << q << " = -2 *( L0(" << maxL0 << ") - maxL(" << maxL <<") )"<< std::endl;
    if ( pL->signalNorm >= 0 && q > 0 ) //this is to catch roundoff error, but it could hide bugs
    {
        return - 2 * ( maxL0 - maxL );  
    }
    else
    {
        //std::cout << "returning q=0 " << pL->signalNorm << " " << maxL0 << " " << maxL << std::endl;
        return 0;
    }
}

//searching for the mu which gives a median significance of 4.38sigma
double my_q0(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    paramList *pL = (paramList *)params;    
    
    for(int i=0; i < pL->source.numFlux; i++)
        pL->source.nuFluxNorm[i] = 1;
        
    pL->detectors[pL->detj].BgNorm = 1;
    pL->signalNorm = gsl_vector_get(v, 0);
    if( pL->signalNorm < 0 )
        return 1e99;
    generateBinnedData( pL, pL->detj, 1, simSeed);
    
    return pow( sqrt(q0(pL)) - 4.28, 2);  //arbitrary function with a minima at 3 sigma 90% of the time
}


double findCoeff3sig(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x,*dx;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = my_q0;
    my_func.params = (void *)pL;

    //start point
    x = gsl_vector_alloc (1);
    dx = gsl_vector_alloc (1);
    
    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/10);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        status = gsl_multimin_fminimizer_iterate (s);
        iter++;
        //std::cout << iter << " " << gsl_vector_get(s->x,0)*pL->C << ", q' = " << s->fval << ", size " << gsl_multimin_fminimizer_size(s) << std::endl;
    }
    while (iter < 30 && s->fval > .0014 && !status); //under 1% error in 4.28 sigma value
        
    double mu = gsl_vector_get(s->x, 0);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    if(iter==30)
    {
        double approxError = sqrt(s->fval)/4.28*100;
        std::cout << "WARNING: non-convergence, sigma - 4.28 = " << s->fval << " " << approxError << "% error" << std::endl;
        if (approxError > 2)
            return 0;
        else
            return mu;
    }
    else
        return mu;
}

void discLimitEvolution(paramList *pL, int detj)
{

    std::cout << "Starting disc. evolution calculations..." << std::endl;

    char filename[100];
    std::ofstream outfile;
    
    sprintf(filename, "%sdiscEvo_xe.dat",pL->root);
    
    std::cout << "writing output to: " << filename << std::endl;    
    outfile.open(filename,std::ios::out);
         
    //determine first guess for mu, need a mu that gives BSM ~ SM/100
    pL->signalNorm=.1;
    
    pL->signalNorm = mu;         
    
    double coup;
    
    while (pL->detectors[detj].exposure < 100)
    {

        mu = findCoeff3sig(pL);

        if (mu==mu) //check for NAN
        {
            //print out result
            std::cout << pL->detectors[detj].exposure << "  " << coup << std::endl;
            outfile   << pL->detectors[detj].exposure << "  " << coup << std::endl;
            
            pL->signalNorm = mu/2;      //update guess
        }
        else
        {
            pL->signalNorm=.001;
            goto makeAguessEvo;
        }
        
        pL->detectors[detj].exposure*=1.2; //increment exposure
        
    }
    outfile.close();

}
