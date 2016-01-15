#ifndef ROBUSTNESS_H
#define ROBUSTNESS_H

#include <gsl/gsl_vector.h>

gsl_vector *EvaluateRobustness(gsl_vector *, int, int, int, gsl_vector*);	

#endif
