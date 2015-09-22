#ifndef BEF_NEW_H
#define BEF_NEW_H_H

#include "structs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>				


double* predOnRes(struct foodweb nicheweb, const double[], double*); 

double* intraguildPred(struct foodweb nicheweb,  const double[], double*);	
												
double* metabolicLoss(struct foodweb nicheweb, const double[], double*);				

double* functionalDiversity(struct foodweb nicheweb, const double[], double[], double[], double*);

int kastenFunktion(double, double, double);		


#endif