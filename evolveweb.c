/*	
	2015-07-02 16:02:21 
	Projekt: Nischennetz mit Migration /Robustness

	Quelltextdatei für komplettes Nischennetz auf Y Patches
 
		Inhalt: 	EvolveNetwork

	Zugehöriger header: evolveweb.h
	
	Diese Funktion berechnet die Populationsdynamik eines gegebenen Nischen-Netzwerkes mit Hilfe eines ODE Solvers.
	Die Dynamik wird durch eine Holling Typ II Funktion charakterisiert. Der Rückgabewert enthält 5*(Rnum+S)*Y + 9 + S Werte:
		
	y	(Rnum+S): Die Größen der Spezies und Ressourcen nach Lösen der DGL
	y0	(Rnum+S): Die Populationsgrößen am Anfang
	ymax(Rnum+S): Die größte Population zwischen t1 und t2
	ymin(Rnum+S): Die kleinste Population zwischen t1 und t2
	yavg(Rnum+S): Die durchschnittliche Population
	Fixp	   6: Fixpunktvariablen
	TL		   S: Trophische Level der Spezies zur Auswertung
*/
#include "structs.h"
#include "holling2.h"
#include "evolveweb.h"
#include "BEF_new_abTL2.h"

#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

gsl_vector* EvolveNetwork(struct foodweb nicheweb, gsl_rng* rng1, const gsl_rng_type* rng1_T)
{	
	struct foodweb *params = &nicheweb; 									// Damit Holling2 auf das foodweb zugreifen kann

	int S 	 	= nicheweb.S;
	int Y 	     	= nicheweb.Y;
	int Rnum 	= nicheweb.Rnum;

	double Rsize = gsl_vector_get(nicheweb.network, (Rnum+S)*(Rnum+S)+Y*Y+2);
	
		
	double *y 	 = (double *)calloc((Rnum+S)*Y, sizeof(double));				// Ergebnis Array für den Lösungsalgorithmus

	int i, j,l   	 = 0;	
	int closezero= 1;			

//-- Ergebnis Variablen-----------------------------------------------------------------------------------------------------------------------------------	
	gsl_vector *result	= gsl_vector_calloc((Rnum+S)*Y*5 + 3 + S +5*Y); 				// y[Simulation], y0, ymax, ymin, yavg, fixp, TL
	gsl_vector *y0		= gsl_vector_calloc((Rnum+S)*Y);						// Startwerte der Populationsgrößen
	gsl_vector *ymax	= gsl_vector_calloc((Rnum+S)*Y);						// Maximalwerte nach t2
	gsl_vector *ymin	= gsl_vector_calloc((Rnum+S)*Y);						// Minimalwerte nach t2
	gsl_vector *yavg	= gsl_vector_calloc((Rnum+S)*Y);						// Durchschnittswert nach t2
// 	gsl_vector *ytest	= gsl_vector_calloc((Rnum+S)*Y);
// 	gsl_vector *ytest2	= gsl_vector_calloc(Y);
	
//--Zufallszahlengenerator für Populationsgrößen----------------------------------------------------------------------------------------------------------

// 	const gsl_rng_type *rng1_T;							// Für zufällige Populationsgröße der Spezies
// 	gsl_rng *rng1;   									
// 	gsl_rng_env_setup();   								
// 	rng1_T = gsl_rng_default;   						// default random number generator (so called mt19937)
// 	//gsl_rng_default_seed = 0;							// default seed for rng
// 	gsl_rng_default_seed=((unsigned)time(NULL));		// random starting seed for rng
// 	rng1 = gsl_rng_alloc(rng1_T);	

//--Erstelle y[] mit Startwerten für die Speziespopulationsgrößen---------------------------------------------------------------------------------------

	  for(j=0; j<Y; j++)								// set initial species size "RANDOM" in each patch
	  {

		for(i=0; i<Rnum; i++)
		 {
			y[j*(Rnum+S)+i] = Rsize;							// Ressourcen Größe pro Patch		
		 }

		for(i=Rnum; i<Rnum+S; i++)
		{
		  if(closezero == 1) y[j*(Rnum+S)+i] = 0.0000001 + (gsl_rng_uniform_pos(rng1)*0.1);
		  else y[j*(Rnum+S)+i] = 0.001     + (gsl_rng_uniform_pos(rng1)*0.1);

			//printf("y0 = %f\n", y[j*(Rnum+S)+i]);
		}
	  }

	printf("Spezies Anfangspopulationen erzeugt\n");

  	gsl_vector_view y_vec = gsl_vector_view_array(y, (Rnum+S)*Y);		
   					       	gsl_vector_memcpy(y0, &y_vec.vector);						//y0 enthält jetzt die so eben bestimmten Startwerte

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
/*ODE: Ordinary Differential Equation  mittlerweile gibt es die odeiv2 systeme -> sollte man vielleicht upgraden
	Hilfe für die alte Version: http://www.inference.phy.cam.ac.uk/pjc51/local/gsl/manual/gsl-ref_25.html	
	Neue Version: https://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html#Ordinary-Differential-Equations
*/  
	const gsl_odeiv_step_type *Solv			= gsl_odeiv_step_rkf45;							// ODE Solver vom Typ RungeKutta 4/5 (siehe Dokumentation)
	      gsl_odeiv_step *s				= gsl_odeiv_step_alloc(Solv,(Rnum+S)*Y);			// Schrittfunktion
	      gsl_odeiv_control *c			= gsl_odeiv_control_y_new(1e-6, 1e-8);			// Kontrollfunktion zur Anpassung der Schrittgröße, um Genuigkeit zu gewährleisten
	      gsl_odeiv_evolve *e			= gsl_odeiv_evolve_alloc((Rnum+S)*Y);			// 
	      gsl_odeiv_system sys			= {Holling2, NULL, (Rnum+S)*Y, params};			// ODE System struct -> dieses wird evolviert

/*--------------------------------------------------------------------------------------------------------------------------------------------------------------- 
gsl_odeiv_system ist ein Datentyp, der ein allgemeines ODE system mit frei wählbaren Parametern enthält. 
Er wird definiert über vier Größen 

(1) eine int funktion f(double t, const double y[], double dydt[], void * params) 
    Sie sollte die Vektorelemente der zu lösenden Funktion in dydt[] speichern für die Argumente (t,y) und die Parameter params enthalten

	-> hier Hol_dynam(double t, const double y[], double ydot[], void *params)

(2) eine funktion "int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params)", die die Jacobi-Matrix enthält. 
    Sie wird von weiter entwickelten Solvern benutzt. Für einfachere Solver muss sie nicht angegeben werden.
	-> hier weggelassen mit NULL

(3) eine Größe size_t dimension, die die Dimension des Gleichungssystems angibt
	-> hier (Rnum+S)*Y, da auf Y Patches jeweils S+Rnum Spezies leben

(4) void * params, einen Pointer auf die Parameter des Systems, diese werden hier über den *params übergeben
    FRAGE: Muss network erst in ein Array geschrieben werden?
*/

//--DGL lösen mit Holling II---------------------------------------------------------------------------------------------------------------------------------------------------
  double t		= 0.0; 				// start time
  double tend1 		= 7800; 
  double tend2		= 8000;				// endtime
  double h		= 1e-5;				// stepwidth

  double countsteps 	= 0;			// Schritte

//  int docheck 		= 0;			

//--Erster Abschnitt bis t1--------------------------------------------------------------------------------------------------------------  
	
  printf("Starte Lösen der Populationsdynamik\n");	
	
  while(t < tend1)					
  { 
	//for(i=0; i<Rnum+S; i++)printf("y[i]%f\n", y[i]);	// Startgrößen
		
	 int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tend1, &h, y);	 			
	/*status = Ergebnis für einen Zeitschritt der DGL Lösung. 
	  HollingII wird mit t, y, params, aufgerufen, die Ergebnisse legt es dann in ydot ab. 
	  Die Werte in ydot werden dann als neue Werte für y verwendet.*/
	
	//	if(status == GSL_SUCCESS)	printf("Status OK\n");

    if(status != GSL_SUCCESS) {
		printf("Fehler beim Lösen der DGL!\n");		
		break;
	   }

    for(i=0; i<(Rnum+S)*Y; i++)
     {
		  if(y[i] < 1.0e-5) 
			 y[i] = 0;			// bei Populationsgrößen kleiner als 10^-5 gilt die Population als ausgestorben
     }
  }

  for(i=0; i < (Rnum+S)*Y; i++)		// Referenzwerte in min und max schreiben = Wert von y nach t = 7800 
  {
    gsl_vector_set(ymax, i, y[i]);
    gsl_vector_set(ymin, i, y[i]);
  }
 
//-- Lösen von t1 bis t2-----------------------------------------------------------------------------------------------------------------------   

//  docheck = 1;			// triggert, dass in HollingII nach Fixpunkt Attraktoren gesucht wird???

// int testf0, testf1, testf2  	= 1;	

	//printf("t=%f\n", t);	

	
  double *metLossAv=(double *) calloc(Y,sizeof(double));
  double *resPredAv=(double *) calloc(Y,sizeof(double));
  double *intraPredAv=(double *) calloc(Y,sizeof(double));
  double *funcDiversityAv=(double *) calloc(Y,sizeof(double));
  double *intraCompetitionAv=(double *) calloc(Y,sizeof(double));
	
  gsl_vector *network 	= nicheweb.network;						// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
  printf("Versuche auf Element %i zuzugreifen",((Rnum+S)*(Rnum+S))+1+(Y*Y)+1+(Rnum+S)+S+S);
  
  gsl_vector_view nvi_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1+(Rnum+S)+S, S);	// Massenvektor
  gsl_vector *nvivec	   = &nvi_vec.vector;
  
  gsl_vector_view fri_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1+(Rnum+S)+S+S, S);	// Massenvektor
  gsl_vector *frivec	   = &fri_vec.vector;
  
  gsl_vector_view fci_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1+(Rnum+S)+S+2*S, S);	// Massenvektor
  gsl_vector *fcivec	   = &fci_vec.vector;
  
//   for(i=0;i<(Rnum+S)*(Rnum+S);i++)
//   {
//     //printf("index %i von network ist: %f\n",i,gsl_vector_get(network,i));
//   }
  
  while(t < tend2)
  {

	//printf("t=%f\n", t);
    countsteps++;

    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tend2, &h, y);		// Hier werden fixp Variablen benutzt

    if(status != GSL_SUCCESS)
      break;

    double metLoss[Y];
    double intraPred[Y];
    double resPred[Y];
    double funcDiversity[Y];
    double intraCompetition[Y];
    for(i=1;i<S+1;i++)
    {
      //printf("y : %f\n",y[i]);
    }
    
    for(i=0; i<(Rnum+S)*Y; i++)
    {
      if(y[i]< 1.0e-5)				// wieder Aussterbe-Kriterium
      y[i]= 0;
				
    }
      
    double intervall[(Rnum+S)*(Y+1)];
    double center[(Rnum+S)*(Y+1)];
    gsl_vector *intervall_gsl=gsl_vector_calloc(S);	
    
    gsl_vector_memcpy(intervall_gsl,nvivec);
    gsl_vector_mul(intervall_gsl,frivec);
    gsl_vector_scale(intervall_gsl,0.5);
    
    int l;
    for(l = 0; l<Y; l++)
    {
      funcDiversity[l] =0;
      for(i = Rnum; i< Rnum+S;i++)
      {
	intervall[i+l*(Rnum+S)] = gsl_vector_get(intervall_gsl,i-Rnum);
	center[i+l*(Rnum+S)] = gsl_vector_get(fcivec,i-Rnum);
      }
    }
     metabolicLoss(nicheweb, y, metLoss);
     intraguildPred(nicheweb, y, intraPred);
     predOnRes(nicheweb, y, resPred);
     functionalDiversity(nicheweb, y, center, intervall, funcDiversity);
     intraspecificCompetition(nicheweb, y,intraCompetition);
    
    //printf("intraguildPred: %f\n",intraPred[0]);
    for(l=0;l<Y;l++)
    {
      //printf("metLoss von patch %i ist %f\n",l,metLoss[l]);
      metLossAv[l]=((metLossAv[l]*(countsteps-1))+metLoss[l])/countsteps;
      resPredAv[l]=((resPredAv[l]*(countsteps-1))+resPred[l])/countsteps;
      intraPredAv[l]=((intraPredAv[l]*(countsteps-1))+intraPred[l])/countsteps;
      funcDiversityAv[l]=((funcDiversityAv[l]*(countsteps-1))+funcDiversity[l])/countsteps;
      intraCompetitionAv[l]=((intraCompetitionAv[l]*(countsteps-1))+intraCompetition[l])/countsteps;
    }
    //printf("resPredAv: %f\n",resPredAv[0]);
    for(i=0; i<(Rnum+S)*Y; i++)
      {
		  if(y[i] > gsl_vector_get(ymax, i))				// Checken ob y größer geworden ist
			 		gsl_vector_set(ymax, i, y[i]); 	

		  if(y[i] < gsl_vector_get(ymin, i))				// Checken ob y kleiner geworden ist
			 		gsl_vector_set(ymin, i, y[i]); 	

		  gsl_vector_set(yavg, i, ((gsl_vector_get(yavg, i)*(countsteps-1)+y[i])/countsteps));
      }
      
	// if(status == GSL_SUCCESS)	printf("Status OK\n");

  	
  	//testf0 	= testf0*fixp0;				Holling verwendet nur fix0, 1, 2 und fixp 0, 1,2
	//testf1	= testf1*fixp1;				Da testf0,1,2 vorher 1 sind stehen in test0,1,2 die Werte von fixp0,1,2 
    //testf2	= testf2*fixp2;
   
  }
  
//   gsl_vector_set_zero(ytest2);
//   gsl_vector_memcpy(ytest,yavg);
//   
//   for(i= 0; i<Y; i++)
//   {
//     int j;
//     for(j = Rnum ; j < (Rnum +S); j++ )
//     {
//       gsl_vector_set(ytest2, i, gsl_vector_get(ytest2,i) + gsl_vector_get(ytest,j+(Rnum+S)*i));
//     }
//   }
//   
//   printf("ytest2 ist %f\n",gsl_vector_get(ytest2,1));
//   printf("mLoss ist %f\n", metLossAv[1]);
  /* --> hier stimmt die Beziehnung von mLoss und biomasse noch */  
  
//--Ergebnis zusammen fassen--------------------------------------------------------------------------------------------------------------------   

  for(i=0; i<(Rnum+S)*Y; i++)
 	gsl_vector_set(result, 0*Y*(Rnum+S)+i, y[i]);						//y[Ende] + y0 + ymax + ymin + yavg + fixp + TL
	 			 
  for(i=0; i<(Rnum+S)*Y; i++)
	gsl_vector_set(result, 1*Y*(Rnum+S)+i, gsl_vector_get(y0, i));	

  for(i=0; i<(Rnum+S)*Y; i++)
	gsl_vector_set(result, 2*Y*(Rnum+S)+i, gsl_vector_get(ymax, i));	

  for(i=0; i<(Rnum+S)*Y; i++)
	gsl_vector_set(result, 3*Y*(Rnum+S)+i, gsl_vector_get(ymin, i));	

  for(i=0; i<(Rnum+S)*Y; i++)
	gsl_vector_set(result, 4*Y*(Rnum+S)+i, gsl_vector_get(yavg, i));

  for(i=0; i < S; i++)
	gsl_vector_set(result, 5*Y*(Rnum+S)+i, gsl_vector_get(nicheweb.network, (Rnum+S)*(Rnum+S)+1+Y*Y+1+(Rnum+S)+i));
	

  
  // Fixpunkte mit übergeben, die sollen in .out datei
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+0, gsl_vector_get((&nicheweb)->fixpunkte, 3));	
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+1, gsl_vector_get((&nicheweb)->fixpunkte, 4));
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+2, gsl_vector_get((&nicheweb)->fixpunkte, 5));
	
	
  for(l=0; l < Y; l++)
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+3+l, metLossAv[l]);
  
  for(l=0; l < Y; l++)
  {
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+3+Y+l, resPredAv[l]);
	//printf("result: %f\n",gsl_vector_get(result,5*Y*(Rnum+S)+S+3+Y));
  }
  
  for(l=0; l < Y; l++)
  {
    gsl_vector_set(result, 5*Y*(Rnum+S)+S+3+2*Y+l, intraPredAv[l]);
    //printf("intraPredAv: %f\n",intraPredAv[l]);
  }
  
  for(l=0; l < Y; l++)
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+3+3*Y+l, funcDiversityAv[l]);
  
  for(l=0; l < Y; l++)
	gsl_vector_set(result, 5*Y*(Rnum+S)+S+3+4*Y+l, intraCompetitionAv[l]);
  
  //printf("funcDiversityAv: %f\n",gsl_vector_get(result,5*Y*(Rnum+S)+S+3+3*Y+0));
	free(y);
//	free(params);

	gsl_vector_free(y0);
	gsl_vector_free(ymax);
	gsl_vector_free(ymin);
	gsl_vector_free(yavg);
// 	gsl_vector_free(nvivec);
// 	gsl_vector_free(frivec);
// 	gsl_vector_free(fcivec);

 	//gsl_rng_free(rng1);

	gsl_odeiv_control_free(c);
  	gsl_odeiv_evolve_free(e);
	
	
  return result;

}

