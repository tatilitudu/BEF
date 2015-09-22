/*	
	2015-07-06 17:32:16 
	Projekt: Nischennetz mit Migration /Robustness

	Quelltextdatei für komplettes Nischennetz auf Y Patches
 
		Inhalt: 	HollingII(double t, const double y[], double ydot[], void *params)

	Zugehöriger header: holling2.h

	Diese Funktion berechnet eine Populationsdynamik mit der Holling Typ II Form. 
	Die Reihenfolge und Datentypen der Paramter sind vorgegeben durch den ODE Solver aus der Funktion EvolveNetwork. 
	Die Ergebnisse der Berechnung liegen in ydot, diese werden vom Solver an das rufende Programm übergeben, sodass Holling2 selbst nichts zurück gibt.
*/
#include "structs.h"
#include "holling2.h"

#include <string.h>						// string modification functions
#include <time.h>						// time functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_sort_vector.h>		// vector sorting functions
#include <gsl/gsl_odeiv.h>				// differential equation solver
#include <gsl/gsl_errno.h>				// errorhandler

int Holling2(double t, const double y[], double ydot[], void *params){
				


	int i, j,l	= 0;						// Hilfsvariablen
	double rowsum	= 0;	
	double colsum	= 0;		  

//-- Struktur zerlegen-------------------------------------------------------------------------------------------------------------------------------


  	struct foodweb *nicheweb = (struct foodweb *)params;			// pointer cast from (void*) to (struct foodweb*)
	
	gsl_vector *network = (nicheweb->network);						// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S

	int S 	 	= nicheweb->S;
	int Y 	 	= nicheweb->Y;
	int Rnum	= nicheweb->Rnum;
	double d  	= nicheweb->d;
	double dij 	= pow(10, d);
	double alpha	= nicheweb->alpha;				// respiration
	double lambda	= nicheweb->lambda;				// ecologic efficiency
	double hand	= nicheweb->hand;				// handling time
	double beta	= nicheweb->beta;				// intraspecific competition
	double aij	= nicheweb->aij;				// attack rate
	
// 	printf("alpha %f\n",alpha);	
// 	printf("lambda %f\n",lambda);
// 	printf("hand %f\n",hand);
// 	printf("beta %f\n",beta);
// 	printf("aij %f\n",aij);
	
	//int len	 = (Rnum+S)*(Rnum+S)+2+Y*Y+(Rnum+S)+S;
	
	gsl_vector_view A_view = gsl_vector_subvector(network, 0, (Rnum+S)*(Rnum+S));						// Fressmatrix A als Vektor
	gsl_matrix_view EA_mat = gsl_matrix_view_vector(&A_view.vector, (Rnum+S), (Rnum+S));				// A als Matrix_view
	gsl_matrix *EAmat	   = &EA_mat.matrix;															// A als Matrix

	gsl_vector_view D_view = gsl_vector_subvector(network, (Rnum+S)*(Rnum+S)+1, Y*Y);					// Migrationsmatrix D als Vektor
	gsl_matrix_view ED_mat = gsl_matrix_view_vector(&D_view.vector, Y, Y);								// D als Matrixview
	gsl_matrix *EDmat	   = &ED_mat.matrix;															// D als Matrix

	gsl_vector_view M_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
	gsl_vector *Mvec	   = &M_vec.vector;

			
// 			printf("\ncheckpoint Holling2 I\n");
// 			printf("\nS = %i\n", S);
// 			printf("\nS + Rnum = %i\n", S+Rnum);
// 
// 			printf("\nSize A_view = %i\n", (int)A_view.vector.size);
// 			printf("\nSize D_view = %i\n", (int)D_view.vector.size);
// 			printf("\nSize M_vec  = %i\n", (int)M_vec.vector.size);


// 			for(i=0; i<(Rnum+S)*Y; i++){
// 				printf("\ny = %f\n", y[i]);
// 				}

// 			for(i=0; i<(Rnum+S)*Y; i++){
// 			printf("\nydot = %f\n", ydot[i]);
// 			}
		

//--zusätzliche Variablen anlegen-------------------------------------------------------------------------------------------------------------

  double ytemp[(Rnum+S)*Y];		 
	for(i=0; i<(Rnum+S)*Y; i++) ytemp[i] = y[i];							// temp array mit Kopie der Startwerte
 	
  for(i=0; i<(Rnum+S)*Y; i++) ydot[i] = 0;									// Ergebnis, in das evolve_apply schreibt
 						
  gsl_vector_view yfddot_vec = gsl_vector_view_array(ydot, (Rnum+S)*Y);		//Notiz: vector_view_array etc. arbeiten auf den original Daten der ihnen zugeordneten Arrays/Vektoren!
  gsl_vector *yfddotvec		 = &yfddot_vec.vector;							// zum einfacheren Rechnen ydot über vector_view_array ansprechen
  
  gsl_vector_view yfd_vec	 = gsl_vector_view_array(ytemp, (Rnum+S)*Y);
  gsl_vector *yfdvec		 = &yfd_vec.vector;								// Startwerte der Populationen

//-- neue Objekte zum Rechnen anlegen--------------------------------------------------------------------------------------------------------

  gsl_matrix *AFgsl	= gsl_matrix_calloc(Rnum+S, Rnum+S);	// matrix of foraging efforts
  gsl_matrix *ADgsl	= gsl_matrix_calloc(Y,Y); 				// matrix of migration efforts
  
  gsl_matrix *Emat	= gsl_matrix_calloc(Rnum+S, Rnum+S);	// gsl objects for calculations of populations 
  gsl_vector *tvec	= gsl_vector_calloc(Rnum+S);
  gsl_vector *rvec	= gsl_vector_calloc(Rnum+S);
  gsl_vector *svec	= gsl_vector_calloc(Rnum+S);
  
  gsl_matrix *Dmat	= gsl_matrix_calloc(Y,Y);				// gsl objects for calculations of migration
  gsl_vector *d1vec	= gsl_vector_calloc(Y);
  gsl_vector *d2vec	= gsl_vector_calloc(Y);
  gsl_vector *d3vec	= gsl_vector_calloc(Y);
  
//	printf("\ncheckpoint Holling2 III\n");

//-- Einzelne Patches lösen------------------------------------------------------------------------------------------------------------    
  for(l=0; l<Y; l++)								// start of patch solving
  {
    gsl_matrix_set_zero(AFgsl);						// Objekte zum Rechnen vor jedem Patch nullen 
    gsl_matrix_set_zero(Emat);
    gsl_vector_set_zero(tvec);
    gsl_vector_set_zero(rvec);
    gsl_vector_set_zero(svec);
    
    gsl_vector_view ydot_vec = gsl_vector_subvector(yfddotvec, (Rnum+S)*l, (Rnum+S));	// enthält ydot von Patch l
    gsl_vector *ydotvec 	 = &ydot_vec.vector;

    gsl_vector_view y_vec	 = gsl_vector_subvector(yfdvec, (Rnum+S)*l, (Rnum+S));		// enthält Startwerte der Population in l
    gsl_vector *yvec 		 = &y_vec.vector;
    
    gsl_matrix_memcpy(AFgsl, EAmat);
    
    for(i=0; i<Rnum+S; i++)
    {
      gsl_vector_view rowA   = gsl_matrix_row(AFgsl,i);
      				  rowsum = gsl_blas_dasum(&rowA.vector);
      if(rowsum !=0 )
      {
		for(j=0; j<Rnum+S; j++)
	    gsl_matrix_set(AFgsl, i, j, (gsl_matrix_get(AFgsl,i,j)/rowsum));				// normiere Beute Afgsl = A(Beutelinks auf 1 normiert) = f(i,j)
      }
    }
    
    gsl_matrix_memcpy(Emat, EAmat);									//  Emat = A
    gsl_matrix_scale(Emat, aij);									//  Emat(i,j) = a(i,j)
    gsl_matrix_mul_elements(Emat, AFgsl);							//  Emat(i,j) = a(i,j)*f(i,j)

    gsl_vector_memcpy(svec, yvec);									// s(i) = y(i)
    gsl_vector_scale(svec, hand);									// s(i) = y(i)*h
    gsl_blas_dgemv(CblasNoTrans, 1, Emat, svec, 0, rvec);			// r(i) = Sum_k h*a(i,k)*f(i,k)*y(k)
    gsl_vector_add_constant(rvec, 1);								// r(i) = 1+Sum_k h*a(i,k)*f(i,k)*y(k)
    	
    gsl_vector_memcpy(tvec, Mvec);									// t(i) = masse(i)^(-0.25)
    gsl_vector_div(tvec, rvec);										// t(i) = masse(i)^(-0.25)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))
    gsl_vector_mul(tvec, yvec);										// t(i) = masse(i)^(-0.25)*y(i)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))

    gsl_blas_dgemv(CblasTrans, 1, Emat, tvec, 0, rvec);				// r(i) = Sum_j a(j,i)*f(j,i)*t(j)
    gsl_vector_mul(rvec, yvec);										// r(i) = Sum_j a(j,i)*f(j,i)*t(j)*y(i) [rvec: Praedation]

    gsl_blas_dgemv(CblasNoTrans, lambda, Emat, yvec, 0, ydotvec);	// ydot(i) = Sum_j lambda*a(i,j)*f(i,j)*y(j)
    gsl_vector_mul(ydotvec, tvec);									// ydot(i) = Sum_j lambda*a(i,j)*f(i,j)*y(j)*t(i)
    
    gsl_vector_memcpy(svec, Mvec);
    gsl_vector_scale(svec, alpha);								// s(i) = alpha*masse^(-0.25) [svec=Respiration bzw. Mortalitaet]

    gsl_vector_memcpy(tvec, Mvec);
    gsl_vector_scale(tvec, beta);								// t(i) = beta*masse^(-0.25)
    gsl_vector_mul(tvec, yvec);									// t(i) = beta*y(i)
    gsl_vector_add(svec, tvec);									// s(i) = alpha*masse^(-0.25)+beta*y(i)
    	
    gsl_vector_mul(svec, yvec);									// s(i) = alpha*masse^(-0.25)*y(i)+beta*y(i)*y(i)
    gsl_vector_add(svec, rvec);									// [svec: Respiration, competition und Praedation]
    
    gsl_vector_sub(ydotvec, svec);								// ydot(i) = Fressen-Respiration-Competition-Praedation
    
    for(i=0; i<Rnum; i++)
      gsl_vector_set(ydotvec, i, 0.0);							// konstante Ressourcen
      
  }// Ende Einzelpatch, Ergebnis steht in ydotvec 

//	printf("\ncheckpoint Holling2 IV\n");
  
//-- Migration lösen---------------------------------------------------------------------------------------------------------    
  for(l = Rnum; l< Rnum+S; l++)								// start of migration solving
  {
    gsl_matrix_set_zero(ADgsl);								// reset gsl objects for every patch
    gsl_matrix_set_zero(Dmat);    
    gsl_vector_set_zero(d1vec);
    gsl_vector_set_zero(d2vec);
    gsl_vector_set_zero(d3vec);

	// Untervektor von yfddot (enthält ydot[]) mit offset l (Rnum...Rnum+S) und Abstand zwischen den Elementen (stride) von Rnum+S.
	// Dies ergibt gerade die Größe einer Spezies in jedem Patch in einem Vektor
    gsl_vector_view dydot_vec = gsl_vector_subvector_with_stride(yfddotvec, l, (Rnum+S), Y);	// ydot[]		
    gsl_vector *dydotvec	  = &dydot_vec.vector;

    gsl_vector_view dy_vec	  = gsl_vector_subvector_with_stride(yfdvec, l, (Rnum+S), Y);			// Startgrößen der Spezies pro Patch
    gsl_vector *dyvec		  = &dy_vec.vector;
          
    gsl_matrix_memcpy(ADgsl, EDmat);		// ADgsl = D
    
    if(nicheweb->M == 1)				// umschalten w: patchwise (Abwanderung aus jedem Patch gleich), sonst linkwise (Abwanderung pro link gleich) 
	   {
		  for(i=0; i<Y; i++)
		   {
				gsl_vector_view colD = gsl_matrix_column(ADgsl, i);					// Spalte i aus Migrationsmatrix
							  colsum = gsl_blas_dasum(&colD.vector);
				if(colsum!=0)
					{
					  for(j=0;j<Y;j++)
					  gsl_matrix_set(ADgsl,j,i,(gsl_matrix_get(ADgsl,j,i)/colsum));		// ADgsl: D mit normierten Links
					}
		    }
	   }

    gsl_matrix_memcpy(Dmat, EDmat);								// Dmat = D
    gsl_matrix_scale(Dmat, dij);								// Dmat(i,j) = d(i,j) (Migrationsstärke)
    gsl_matrix_mul_elements(Dmat, ADgsl);						// Dmat(i,j) = d(i,j)*xi(i,j)   (skalierte und normierte Migrationsmatrix)
     
    gsl_vector_set_all(d1vec, 1/gsl_vector_get(Mvec, l));		// d1(i)= m(l)^0.25
    gsl_vector_mul(d1vec, dyvec);								// d1(i)= m(l)^0.25*y(i)
    gsl_blas_dgemv(CblasNoTrans, 1, Dmat, d1vec, 0, d2vec);		// d2(i)= Sum_j d(i,j)*xi(i,j)*m(l)^0.25*y(j)
    
    gsl_vector_set_all(d1vec, 1);								// d1(i)= 1
    gsl_blas_dgemv(CblasTrans, 1, Dmat, d1vec, 0, d3vec);		// d3(i)= Sum_j d(i,j)*xi(i,j)
    gsl_vector_scale(d3vec, 1/gsl_vector_get(Mvec,l));			// d3(i)= Sum_j d(i,j)*xi(i,j)*m(l)^0.25
    gsl_vector_mul(d3vec, dyvec);								// d3(i)= Sum_j d(i,j)*xi(i,j)*m(l)^0.25*y(i)
    
    gsl_vector_add(dydotvec, d2vec);				// 
    gsl_vector_sub(dydotvec, d3vec);				// Ergebnis in dydotvec (also ydot[]) = Sum_j d(i,j)*xi(i,j)*m(l)^0.25*y(j) - Sum_j d(i,j)*xi(i,j)*m(l)^0.25*y(i) 
  }	// Patch i gewinnt das was aus allen j Patches zuwandert und verliert was von i auswandert


	//printf("\ncheckpoint Holling2 V\n");

	/*
	for(i=0; i<(Rnum+S)*Y; i++){
		printf("\ny = %f\tydot=%f\n", y[i], ydot[i]);
		}
    */
//--check for fixed point attractor-----------------------------------------------------------------------------------
	
	if(t>7800){

		gsl_vector_set(nicheweb->fixpunkte, 0, 0);	
		gsl_vector_set(nicheweb->fixpunkte, 1, 0);
		gsl_vector_set(nicheweb->fixpunkte, 2, 0);		 

		int fix0 = (int)gsl_vector_get(nicheweb->fixpunkte, 0);
		int fix1 = (int)gsl_vector_get(nicheweb->fixpunkte, 1);
		int fix2 = (int)gsl_vector_get(nicheweb->fixpunkte, 2);


//	printf("t = %f\n", t);
	
		for(i=0; i<(Rnum+S)*Y; i++)
		  {
			  if(y[i] <= 0)
			  {
				fix0++;
				fix1++;
				fix2++;
			  }
			  else 
			  {
				if((ydot[i]/y[i]<0.0001) || (ydot[i]<0.0001)) fix0++;
				if(ydot[i]/y[i]<0.0001) fix1++;
				if(ydot[i]<0.0001) fix2++;
			  }
		  }

    if(fix0==(Rnum+S)*Y) gsl_vector_set(nicheweb->fixpunkte, 3, 1);
    if(fix1==(Rnum+S)*Y) gsl_vector_set(nicheweb->fixpunkte, 4, 1);
    if(fix2==(Rnum+S)*Y) gsl_vector_set(nicheweb->fixpunkte, 5, 1);
  }

//--Speicher leeren----------------------------------------------------------------------------------------------------- 

  gsl_matrix_free(Emat);  
  gsl_matrix_free(Dmat);  
  gsl_matrix_free(AFgsl);  
  gsl_matrix_free(ADgsl);
  
  gsl_vector_free(tvec);
  gsl_vector_free(rvec);
  gsl_vector_free(svec);
  gsl_vector_free(d1vec);
  gsl_vector_free(d2vec);
  gsl_vector_free(d3vec);
  
//	printf("\nCheckpoint Holling2 VI\n");

  return GSL_SUCCESS;

}


