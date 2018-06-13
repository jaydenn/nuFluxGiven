//definition of detector struct
#ifndef STDIO_H
    #include <stdio.h>
#endif
#ifndef GSL_SPLINE_H
    #include <gsl/gsl_spline.h>
#endif
#define DETECTORSTRUCT_H
#define MAXBINS 1
#define INTERP_POINTS 10000

struct detector {
    char name[20];
    int sourcej; //index of source for this detector
    double exposure;
    double AM;
    int nIso;
    int isoZ[10];
    int isoA[10];
    double **ionization;
    double isoFrac[10];
    double isoSZ[10];
    double isoSN[10];
    double isoJN[10];
    double ErL;
    double ErU;
    
    int bg;
    int eff;
    int res;
    
    double *binnedData;     //array of i bins with binnedData[i] number of events per bin
    int nbins;
    double *binW;
    double *unbinnedData;   //array of i events which occured at energy unbinnedData[i]
    double nEvents;
    
    double BgNorm, BgUn, krBgNorm, rnBgNorm, xeBgNorm;  //norm factor for background and fractional uncertainty
    
    gsl_spline *background;
    gsl_interp_accel *accelBg;
    gsl_spline *krBackground;
    gsl_interp_accel *accelkrBg;
    gsl_spline *rnBackground;
    gsl_interp_accel *accelrnBg;
    gsl_spline *xeBackground;
    gsl_interp_accel *accelxeBg;

    gsl_spline *NR_O;
    gsl_interp_accel *accelO;   
    
    gsl_spline *NR_N;
    gsl_interp_accel *accelN;   
    
    gsl_spline *NR_PP;
    gsl_interp_accel *accelPP;
    
    gsl_spline *NR_PEP;
    gsl_interp_accel *accelPEP;
    
    gsl_spline *NR_BE;
    gsl_interp_accel *accelBE;
     
    gsl_spline *NR_B;
    gsl_interp_accel *accelB; 
    
    gsl_spline *NR_F;
    gsl_interp_accel *accelF;
    
    void printDetSpecs()
    {
        printf(" %s\n",name);
        
        for (int i=0; i<nIso; i++)
        {
            printf("    Isotope %d - %2.1f%%\n",i+1,isoFrac[i]*100);
            printf("     A  = %d\n",isoA[i]);    
        }
        printf("   source: %d",sourcej);
        printf("   %.1f < Er < %.1f\n",ErL,ErU);
        printf("   bg  = %d\n",bg);
        printf("   eff = %d\n",eff);
        printf("   res = %d\n",res);
    }
    
    detector()
    {
        nIso=0; AM=-1; isoA[0]=-1; isoFrac[0]=-1; ErL=0; ErU=-1; bg=-1; eff=-1; res=-1, nEvents=0;
        BgNorm = 1; BgUn = -1e-99;
        
        ionization = new double*[10]();
        for(int i=0;i<10;i++)
            ionization[i] = new double [100]();

        background = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelBg = gsl_interp_accel_alloc();
        
        NR_O = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelO = gsl_interp_accel_alloc(); 
        
        NR_N = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelN = gsl_interp_accel_alloc();   
        
        NR_PP = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelPP = gsl_interp_accel_alloc();
        
        NR_PEP = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelPEP = gsl_interp_accel_alloc();
    
        NR_BE = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelBE = gsl_interp_accel_alloc();
        
        NR_B = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelB = gsl_interp_accel_alloc();
        
        NR_F = gsl_spline_alloc(gsl_interp_linear,INTERP_POINTS);
        accelF = gsl_interp_accel_alloc();
    
    
    }
    
    
};

