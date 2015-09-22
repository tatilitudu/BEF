#ifndef HOLLING2_H
#define HOLLING2_H

#include <string.h>						// string modification functions
#include <time.h>						// time functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <gsl/gsl_sort_vector.h>		// vector sorting functions
#include <gsl/gsl_odeiv.h>				// differential equation solver
#include <gsl/gsl_errno.h>				// errorhandler
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

int Holling2(double, const double[], double[], void *);

#endif
