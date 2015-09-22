#ifndef STRUCTS_H
#define STRUCTS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
	Die foodweb Struktur enthält die Konsolenparameter S, B, Rnum, Y, T, M, d, x sowie einen Netzwerk Vektor für die Ergebnisse der Nischennetz-Berechnung.
*/
struct foodweb{

	gsl_vector* network;
	gsl_vector* fixpunkte;				// fix0,1,2, fixp0,1,2, testf0,1,2

	int S;
	int B;
	int Rnum;

	int Y;
	int T;
	
	double d;
	double x;

	double alpha;
	double hand;
	double beta;
	double lambda;
	double aij;

	int M;
	
};


struct resource{

	double size;
	double growth;

};



#endif
