#ifndef MIGNICHEWEB_H
#define MIGNICHEWEB_H

#include "structs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>				


gsl_vector* SetNicheNetwork(struct foodweb, struct resource);

gsl_matrix* SetNicheValues(struct foodweb, double, gsl_rng*, const gsl_rng_type*);	
												
gsl_matrix* SetFeedingMatrix(struct foodweb, gsl_matrix*, double, double);				

gsl_matrix* SetMasses(struct foodweb, gsl_matrix*, gsl_matrix*, double);

gsl_vector* LinkElements(struct foodweb, gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, double, int);		

int CountLinks(gsl_matrix*, int);	

#endif
