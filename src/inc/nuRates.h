#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double nuRatesInit(paramList *pL, int ndet);

double intNuRate(double ERmin, double ERmax, paramList *pL, int detj);

double diffNuRate(double ER, paramList *pL, int detj);
