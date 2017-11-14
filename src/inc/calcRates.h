int calcRates(paramList *pList);

void rateInit( paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline);
