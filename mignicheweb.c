/*	
	2015-07-02 16:02:21 
	Projekt: Nischennetz mit Migration /Robustness

	Quelltextdatei für komplettes Nischennetz auf Y Patches
 
		Inhalt: 	SetNicheNetwork(int, int, int, int, int)
					LinkElements(gsl_matrix*, gsl_matrix*, gsl_matrix*, int, int, int, double)
					CountLinks(gsl_matrix*, int)

	Zugehöriger header: mignicheweb.h
*/

#include "structs.h"				// foodweb Struktur
#include "mignicheweb.h"
#include "topology.h"

#include <time.h>
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_sort_vector.h>			// vector sorting functions

/*
Diese Funktion erstellt ein Nahrungsnetz nach dem Nischenmodell für S Spezies auf Y Lebensräumen, die über die Kopplungskonstante d verbunden sind. T gibt die Topologie an (siehe SetTopology). Rnum ist die Anzahl der Ressource(n) und Rsize die Startgröße der Ressource(n). C beschreibt die Konnektivität des Netzes, wobei CRange die erlaubte max. Abweichung von C angibt.
*/
gsl_vector *SetNicheNetwork(struct foodweb nicheweb, struct resource res, gsl_rng* rng1, const gsl_rng_type* rng1_T){

int len			= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S+3*(nicheweb.S));	// Länge des Rückabewerts

double C		= 0.15;		// vorgegebenes C
double CRange	= 0.01;

double Rsize	= res.size;

  gsl_vector *result		= gsl_vector_calloc(len);
  gsl_matrix *NV		= gsl_matrix_calloc(3, nicheweb.S);
  gsl_matrix *A			= gsl_matrix_calloc((nicheweb.Rnum+nicheweb.S),(nicheweb.S+nicheweb.Rnum));
  gsl_matrix *mas		= gsl_matrix_calloc(2, (nicheweb.Rnum+nicheweb.S));

  printf("Erzeuge Nischennetz mit %i Spezies (davon gewünschte Anzahl basaler Spezies %i) auf %i Patches. Topologie-Marker %i.\n", nicheweb.S, nicheweb.B, nicheweb.Y, nicheweb.T);

//--Zufallszahlengenerator initialisieren--------------------------------------------------------------------------------

// 	const gsl_rng_type *rng1_T;											// ****
// 	gsl_rng *rng1;   													// initialize random number generator
// 	gsl_rng_env_setup();   												// ermöglicht Konsolenparameter
// 	rng1_T = gsl_rng_default;   										// default random number generator (so called mt19937)
// 	gsl_rng_default_seed = 0;											// default seed for rng
// 	//gsl_rng_default_seed = ((unsigned)time(NULL));	// random starting seed for rng
// 	//gsl_rng_default_seed = (rdtsc());
// 	rng1 = gsl_rng_alloc(rng1_T);


//--Alle benötigten Daten ausrechnen--------------------------------------------------------------------------------------

int flag = 1;
int i 	 = 0;

while(flag == 1)
{

	flag = 0; 

	NV 	= SetNicheValues(nicheweb, C, rng1, rng1_T);							// Nischenwerte		NOTIZ: In den gsl Objekten liegen floats!
	A	= SetFeedingMatrix(nicheweb, NV, C, CRange);							// interaction matrix 

	int links = CountLinks(A, (nicheweb.Rnum+nicheweb.S));	

	//printf("Prüfe Linkanzahl\n");
	double CON = (double)links/(((double)nicheweb.S)*((double)nicheweb.S)-1);

	if(CON < (C-CRange) || CON > (C+CRange))
	{
		flag = 1;
		continue;
	}

	int bcount = 0;
	for(i = nicheweb.Rnum; i< (nicheweb.S + nicheweb.Rnum); i++)
	{
		if(gsl_matrix_get(A, i, 0)==1)
		bcount++;			
	}

	if(bcount!=nicheweb.B)
	{			
		flag = 1;
		continue;
	}	

	mas	= SetMasses(nicheweb, NV, A, Rsize);							// Massen + TL

	int TL0 = 0;
	for(i = 0; i< nicheweb.S; i++){
		if(gsl_matrix_get(mas, 1, i+nicheweb.Rnum) == 0)
		TL0++;			
	}

	if(TL0!=0)
	{			
		flag = 1;
		continue;
	}	

}

	for(i=0; i<nicheweb.S; i++) printf("Nischenwerte: %f\n", gsl_matrix_get(NV, 0, i));	
	
	gsl_matrix *D    = SetTopology(nicheweb.Y, nicheweb.T);								// migration matrix

 	result = LinkElements(nicheweb, NV, A, mas, D, Rsize, len);	
		
 		 printf("\nNetzwerk erfolgreich erzeugt!\n");

  gsl_matrix_free(D);
  gsl_matrix_free(NV);
  gsl_matrix_free(A);
  gsl_matrix_free(mas);

	return result;
}

//######################################################################################################################################
/*
Diese Funktion berechnet Nischenwerte für S Spezies. Dabei wird zunächst jeder Spezies eine Zufallszahl zugeordnet, der Nischenwert. 
Dieser ist gleichverteilt auf ]0,1[. Ausgehend davon wird dann ein Fresszentrum und ein Fressbereich bestimmt. 
Der Fressbereich wird mit einer Beta-Verteilung erwürfelt.
Eine Spezies kann eine andere Spezies fressen, wenn der Nischenwert der Beute im Fressbereich des Räubers liegt.
Rückgabewert: 3xS Matrix mit [0][S]: Nischenwert, [1][S]: Fressbereich, [2][S]: Fresszentrum.
*/
gsl_matrix *SetNicheValues(struct foodweb nicheweb, double C, gsl_rng* rng1, const gsl_rng_type* rng_T){

	int S = nicheweb.S;
	//printf("\nStarte Berechnung der Nischenwerte für %i Spezies\n", S);

	gsl_matrix *NV	= gsl_matrix_calloc(3, nicheweb.S);	
 	gsl_vector *nv 	= gsl_vector_calloc(S); 
  
	//printf("nischenwert allokation");

 	double disbeta = (1-2*C)/(2*C);			// Für den Fressbereich (Beta-Verteilung)
 	int i = 0;


//--Nischenwerte ausrechnen------------------------------------------------------------------------------------------------

	for(i= 0; i<S; i++) 
	  gsl_vector_set(nv, i, gsl_rng_uniform_pos(rng1));				// Nischenwerte gleichverteilt auf ]0,1[
	
	gsl_sort_vector(nv);											// Sortieren für Massenberechnung später

	for(i = 0; i < S; i++)
		{	
			double nvi 	= gsl_vector_get(nv, i);
			double fri 	= gsl_ran_beta(rng1, 1, disbeta);
			double rand	= gsl_rng_uniform_pos(rng1);
			double fci 	= nvi*fri*rand/2 + nvi*(1-rand);				// Zufälliges Fresszentrum in [nv(i)*fr(i)/2, nv(i)] 

		    gsl_matrix_set(NV, 0, i, nvi);		
		    gsl_matrix_set(NV, 1, i, fri);
		    gsl_matrix_set(NV, 2, i, fci);
			
 		}

//--Zuweisung----------------------------------------------------------------------------------------------------
		
  free(nv);
	
	return NV;
 
}//end SetNicheValues

//############################################################################################################################################################
/*
Diese Funktion verknüpft die Nischenwertdaten einzelner Spezies zu einer Fressmatrix. Dabei ist der Zeilenindex der Fressmatrix den Räubern zugeordnet und die Spalte der Beute. 
Ein Null-Eintrag bedeutet, dass keine Fressbeziehung besteht. Eine 1 bedeutet, dass die Zeilenspezies die Spaltenspezies fressen kann.

Die Konnektivität (Anzahl vorhandener Links/ Anzahl möglicher Links) sollte innerhalb von C+-CRange liegen, ansonsten wird ein Hinweis angegeben, dass C zu groß/klein ist.

Für Spezies, die keine anderen Spezies fressen können, wird ein Link zur Ressource eingerichtet. Die Anzahl der basalen Spezies wird zum Programmstart vorgegeben (B).
Gibt es eine abweichende Anzahl von basalen Spezies wird dies dem Benutzer mitgeteilt, damit Anpassungen gemacht werden können.
*/
gsl_matrix *SetFeedingMatrix(struct foodweb nicheweb, gsl_matrix* NV, double C, double CRange){

int S 	 = nicheweb.S;
int Rnum = nicheweb.Rnum;

int i, j 	= 0;
int k		= 0;

double nvi 	= 0;
double fri 	= 0;
double fci 	= 0;
double nvj 	= 0;

 gsl_matrix* A = gsl_matrix_calloc((nicheweb.Rnum+nicheweb.S), (nicheweb.Rnum+nicheweb.S));

	//printf("Starte Berechnung der Fressmatrix...\n");
			
	for(i = 0; i< S; i++)
	{
	  nvi	= gsl_matrix_get(NV, 0, i);			// NV = 3xS Matrix! -> 0,1,2 x 0,...S-1 Felderindizes
	  fri	= gsl_matrix_get(NV, 1, i);
	  fci	= gsl_matrix_get(NV, 2, i);

	  for(j	= 0; j< S; j++)
	  	{
			nvj	= gsl_matrix_get(NV, 0, j);

			if( (fci-nvi*fri/2 < nvj) && (nvj < fci+nvi*fri/2) )		// falls nvj im Fressbereich von i liegt setze Link bei A[i+Rnum][j+Rnum]
			  gsl_matrix_set(A, i+Rnum, j+Rnum, 1);							
		}

//--Kannibalismus vermeiden----------------------------------------------------------------------------------------------------------	
  
	  gsl_matrix_set(A, i+Rnum, i+Rnum, 0);											
	 // printf("Kannibalismus vermeiden: A[%i][%i] gesetzt auf %f\n", i+Rnum, i+Rnum, gsl_matrix_get(A, i+Rnum, i+Rnum));	

//--Basale Spezies finden-------------------------------------------------------------------------------------------------------------
	  
	  k = 0;
	  gsl_vector_view tempp = gsl_matrix_row(A, i+Rnum);	// Falls eine Spezies keine Beute hat (Zeile hat keinen Eintrag != 0) wird sie zur Ressource verlinkt
	  k =  gsl_vector_isnull(&tempp.vector);				// 1: Vektor ist 0-Vektor, 0: Vektor ist nicht der 0-Vektor

	  if(k==1)  
		{ 
			gsl_matrix_set(A, i+Rnum, 0, 1);			// A[i][0] gibt an wer die Ressource fressen darf (Ressource darf sich nicht selbst fressen)
		}
	}
	
	return A;

}//end SetFeedingMatrix


//#################################################################################################################################################
/*
	Diese Funktion ordnet an Hand einer Fressmatrix von S Spezies diesen Spezies trophische Level und Massen zu.	
	Sie gibt eine 2x(Rnum+S) große Matrix zurück. 
	Die erste Zeile enthält die Masse der Rnum+S Spezies. 
	Dabei haben basale Spezies die Masse 0 und nicht basale Spezies erhalten ihren Nischenwert als Masse.
	Die trophischen Level der Spezies werden in der zweiten Zeile der Matrix gespeichert. Basale Spezies haben TL 0.
*/
gsl_matrix *SetMasses(struct foodweb nicheweb, gsl_matrix* NV, gsl_matrix* A, double Rsize){		//Ohne x für Allometrie

	//printf("SetMasses");

int S 			= nicheweb.S;
int Rnum		= nicheweb.Rnum;

gsl_matrix* mas = gsl_matrix_calloc(2, Rnum + S);					// nullte Zeile Masse, erste Zeile trophisches Level	
	 
//double x		= nicheweb.x;

int check		= 0; 
int i, j 		= 0;
int tlgesucht 		= 1;				// TL 0 darf es nicht geben, sonst keine Beute für diese Spezies

//--TL = 0: Ressource--------------------------------------------------------------------------------------------------------

for(i=0; i< Rnum; i++) gsl_matrix_set(mas, 1, i, 0);	// Ressource hat TL = 0;

//--TL = 1 Spezies bestimmen (Link zur Ressource)-------------------------------------------------------------------------------

  for(i=Rnum; i < S+Rnum; i++)
	{
	  for(j=0; j<Rnum; j++)
		{
		  if(gsl_matrix_get(A, i, j)!=0)
		  {
			gsl_matrix_set(mas, 1, i, 1);					
			//printf("Spezies %i hat trophisches Level 1\n", i-Rnum);
		  }
		}
	}

//--TL>1 Spezies bestimmen--------------------------------------------------------------------------------	
  
  tlgesucht = 2;

  while(tlgesucht <= S)
   {
	 check = 0;
	
		 for(i=Rnum; i< (S+Rnum); i++)
		  {

			if(gsl_matrix_get(mas, 1, i) == 0)			//noch kein TL zugewiesen -> suche mögliche Beuten
			  {
				for(j=Rnum; j< (S+Rnum); j++)
				{
					//printf("check=%i\n", check);
					//printf("i=%i,\tj=%i\n", i, j);

				  //Spezies hat trophisches Level, das genau 1 kleiner ist und ist auch Beute laut A

				  if((gsl_matrix_get(mas, 1, j) == (tlgesucht-1)) && (gsl_matrix_get(A, i-Rnum, j-Rnum)!=0))		
					  
					  {
						gsl_matrix_set(mas, 1, i, tlgesucht);
						//printf("Spezies %i hat trophisches Level %i\n", i-Rnum, tlgesucht);
						check++;  
					  }
				}
			  }
		   } 

	if(check == 0) 				// Keine Zuweisung in der letzten Runde, aufhören
		{
		  tlgesucht = S+1;
		}						
	else tlgesucht++;
	   
   }							// Frage mit TL=5 als Maximalwert?

	int TL0 = 0;
	int TL1 = 0;
	int TL2 = 0;
	int TL3 = 0;
	int TL4 = 0;
	int TL5 = 0;
	
	for(i=Rnum; i<S+Rnum; i++)
		{
			int TL = (int)gsl_matrix_get(mas, 1, i);
			
			switch(TL){
						case 0:	 TL0++;
								 break;	
						
						case 1:	 TL1++;
								 break;	

						case 2:	 TL2++;
								 break;	

						case 3:	 TL3++;
								 break;	
				
						case 4:	 TL4++;
								 break;	
					
						default: TL5++;
					  }

		}	

	printf("Trophische Level-Verteilung: Ressourcen: %i TL1: %i, TL2: %i, TL3: %i, TL4: %i, TL>4: %i\n", TL0, TL1, TL2, TL3, TL4, TL5);

//--Massen eintragen (0 oder Nischenwert, Rsize bei Ressource)----------------------------------------------------------------------------------

  	gsl_matrix_set(mas, 0, 0, Rsize);				// Nur für eine Ressource gültig! mas[0][0] Größe der Ressource

  for(i = Rnum; i< Rnum+S; i++)
	{ 
	  /*if(nicheweb.x == 0.0){
		  
		   if(gsl_matrix_get(mas, 1, i) == 1)		gsl_matrix_set(mas, 0, i, 0);								// Basale Sp. Masse = 0
		   else if(gsl_matrix_get(mas, 1, i)  > 1) 	gsl_matrix_set(mas, 0, i, gsl_matrix_get(NV, 0, i-Rnum));	// Sonst Masse = Nischenwert 
		  }
	  else{*/
			gsl_matrix_set(mas, 0, i, pow(10, -0.25*(nicheweb.x*4*gsl_matrix_get(NV, 0, i-Rnum))));			// Allometrie
		  //}

	}

	return mas;
}
//end SetMasses


//#####################################################################################################################################################
/*
	LinkElements verbindet die zuvor berechneten Daten eines Nahrungsnetzes in einen Vektor, damit die Populationsdynamik damit berechnet werden kann.
	Die Größe des Vektors beträgt ((Rnum+S)²+ 1 + Y² + 1 + (Rnum+S) + S ). Beachte dabei, dass der Index von 0 bis (Rnum+S)²+...+ S - 1 läuft! 
	Der Rückgabewert enthält:

	(Rnum*S)²: Zeilen von A aneinander gereiht
		  + 1: Links in A
		  Y² : Zeilen von D aneinander gereiht
		  + 1: Links in D
	   Rnum+S: Massen aller Spezies inkl. Ressourcen
	        S: Trophische Level aller Spezies
*/
gsl_vector* LinkElements(struct foodweb nicheweb, gsl_matrix* NV, gsl_matrix* A, gsl_matrix* mas, gsl_matrix* D, double Rsize, int len){

int S 	 = nicheweb.S;
int Y 	 = nicheweb.Y;
int Rnum = nicheweb.Rnum;

int i, j 	= 0;
int index 	= 0;										// Läuft die Elemente von result ab

gsl_vector* result = gsl_vector_calloc(len); 
	
//--Inhalt an die richtigen Stellen schreiben--------------------------------------------------------------------
//--A und Links in A---------------------------------------------------------------------------------------------

	for(i=0; i< (S+Rnum); i++)	
	  {
		for(j=0; j< (S+Rnum); j++)
	 	{
			gsl_vector_set(result, i*(Rnum+S)+j, gsl_matrix_get(A, i, j));	// i*Res+S Zeilen = (Res+S)² Elemente aus A
			//printf("Index = %i, result %i gesetzt auf %f\n", index, i*(Rnum+S)+j, gsl_vector_get(result, i*(Rnum+S)+j));
			index++;
	   	}
		//printf("\n");
	  }

	printf("Fressmatrix im Netzwerk\n");

	gsl_vector_set(result, index, CountLinks(A, Rnum+S));
	//printf("Index = %i, result %i gesetzt auf %f\n\n", index, index, gsl_vector_get(result, index+1));
	index++;

//--D und Links in D-----------------------------------------------------------------------------------------------

	for(i=0; i< Y; i++)			
	  {
		for(j=0; j< Y; j++)
	 	{
			gsl_vector_set(result, index, gsl_matrix_get(D, i, j));	// i*(Res+S) + Y Zeilen = Y² Elemente aus A
			//printf("Index = %i, result %i gesetzt auf %f\n", index, index, gsl_vector_get(result, index));
			index++;
		}
	  }

	printf("\nMigrationsmatrix im Netzwerk\n");

	gsl_vector_set(result, index, CountLinks(D, Y));		//Nochmal checken ob decmigcount wirklich Links in D ist -> JA
	//printf("result %i gesetzt auf %f\n\n", index, gsl_vector_get(result, index));
	index++;

//--Massen auch von Ressource---------------------------------------------------------------------------------------

	gsl_vector_set(result, index, Rsize);
	index++;
	
	for(i=1; i< Rnum+S; i++)	
	  {
	 	gsl_vector_set(result, index, gsl_matrix_get(mas, 0, i));						// i*(Res+S) + Y Zeilen = Y² Elemente aus A
		// printf("result %i gesetzt auf %f\n", index, gsl_vector_get(result, index));
		index++;
	  }

	printf("\nMassen im Netzwerk\n");
	
//--Trophische Level nur von S Spezies-------------------------------------------------------------------------------

	printf("Trophische Level:\n");

	for(i=1; i<= S; i++)
	  {
		gsl_vector_set(result, index, gsl_matrix_get(mas, 1, i));
		printf("%3.1f\t", gsl_vector_get(result, index));
		index++;
	  }

//-- Nischenwerte nur von S Spezies-------------------------------------------------------------------------------



	for(i=0; i< 3; i++)
	  {
	    for(j=0; j<S; j++)
	    {
		gsl_vector_set(result, index, gsl_matrix_get(NV, i, j));
		index++;
	    }
	  }
	  
	printf("Nischenwerte im Netzwerk\n"); 
	printf("\nNetzwerkkomponenten zusammengesetzt. Insgesamt %i Elemente\n\n", index);
	
return result;
 
}// end LinkElements


//#########################################################################################################################
/*
	Zählt die Einträge einer quadratischen Matrix, die nicht 0 sind.
*/
int CountLinks(gsl_matrix* A, int Dim){

int i, j = 0;
int linkCount = 0;


	for(i = 0; i < Dim; i++)
	{
		for(j = 0; j < Dim; j++)
		{
			if(gsl_matrix_get(A, i, j) != 0 )
			linkCount++;
		}
	}

return linkCount;

}// end CountLinks




