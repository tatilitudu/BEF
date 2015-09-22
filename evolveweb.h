#ifndef EVOLVEWEB_H
#define EVOLVEWEB_H

#include "structs.h"
#include "holling2.h"

#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector* EvolveNetwork(struct foodweb);

#endif
