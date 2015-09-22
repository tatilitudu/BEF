/* 2015-07-21 14:50:40 

	Projekt: Nischennetz mit Migration mit Robustness Analyse

	Quelldatei: main.c zu Nischennetz_mit_Migration_1.0
	Programmiersprache: C
	Autoren: S. J. Plitzko, M.Hamm
	Version 1.0

###############################################################################################################################
Beschreibung:
Dieses Programm simuliert ein über mehrere Patches ausgedehntes Nahrungsnetz mit dem Nischenmodell. Die Populationsdynamik wird mit einer Holling Typ II Funktion modelliert. Beim Aufrufen können über die Konsole folgende Parameter bestimmt werden:

Anzahl Spezies S
Anzahl basale Spezies B
Allometrie Koeffizient x
Migrationsstärke d
Topologie des Lebensrausm (siehe Set_Topology())
Anzahl der statistischen Wiederholungen L
Anzahl der Patches Y
Migration ja oder nein "M"

Kompilieren mit: gcc -o V1_15_07_13 main.c getargs.c topology.c mignicheweb.c evolveweb.c holling2.c robustness.c -lm -lgsl -lgslcblas -Wall

Eingabe in der Konsole: ./NAME -S X -B X -T X -L X -Y X -d X.X -x X.X -M "X" -R X.X	mit X bzw. X.X numerischer Wert der Variablen (Reihenfolge egal)

Das Ergebnis der Simulation wird in der Datei xxx.out gespeichert. Dabei sind die Daten wie folgt angeordnet (alles eine Zeile): 
S		B		M		x		Y		dpow	T
Rob		Perlok	Perges
Si_ges	Si_TL1	Si_TL2	Si_TL3	Si_TL4	Si_TL>4		Spezies initial
Sf_ges	Sf_TL1	Sf_TL2	Sf_TL3	Sf_TL4	Sf_TL>4		Spezies final
Bi_ges	Bi_TL1	Bi_TL2	Bi_TL3	Bi_TL4	Bi_TL>4		Biomass initial
Bf_ges	Bf_TL1	Bf_TL2	Bf_TL3	Bf_TL4	Bf_TL>4		Biomass final
Sh_ges	Sh_TL1	Sh_TL2	Sh_TL3	Sh_TL4	Sh_TL>4		Spezies Hub
Bh_ges	Bh_TL1	Bh_TL2	Bh_TL3	Bh_TL4	Bh_TL>4		Biomass Hub
Ss_ges	Ss_TL1	Ss_TL2	Ss_TL3	Ss_TL4	Ss_TL>4		
Bs_ges	Bs_TL1	Bs_TL2	Bs_TL3	Bs_TL4	Bs_TL>4
1mit2	2mit3	3mit1
Fixp0	Fixp1	Fixp2	Fixp3	Fixp4	Fixp5	Fixp6	Fixp7

###############################################################################################################################
*/

//--Standardbibliotheken--------------------------------------------------------------
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

//--Selbstdefinierte Funktionen------------------------------------------------------
#include "structs.h"					// foodweb Struktur

#include "getargs.h"					// getArgs
#include "topology.h"					// SetTopology
#include "mignicheweb.h"				// SetNicheNetwork, LinkElements, CountLinks
#include "evolveweb.h"					// DGL lösen
#include "holling2.h"					// DGL lösen			
#include "robustness.h"					// Robustness Analyse

//--Verzeichnis für Ergebnisse-------------------------------------------------------
#define ORT "/home/tatjana/Arbeitsfläche/MichaelasProgramm/mitBEF/ErsteVersuche/"
//++START++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char** argv)
{	

//--Foodweb Struktur mit Standardwerten aufstellen------------------------------------------------------------------------------------------

	gsl_vector* fixpunkte	= gsl_vector_calloc(9);

	struct foodweb nicheweb	= {NULL, fixpunkte, 18, 3, 1, 5, 0, -7., 0.0, 0.6, 0.35, 0.5, 0.65, 6.0, 0};		// Reihenfolge: network, fxpkt, S, B, Rnum, Y, T, d, x, alpha, hand, beta, lambda, aij, M
	
	struct resource res	= {500.0, 0.0};							// Resource: Größe, Wachstum
	
	
//--Konsoleneingabe-------------------------------------------------------------------------------------------------------------------------
	
	int L = 5;	// Statistik
	int i = 0;	// Counter
	
	int checksum = getArgs(argc, argv, &(nicheweb.S), &(nicheweb.B), &(nicheweb.T), &(nicheweb.d), &L, &(nicheweb.Y), &(nicheweb.x), &(nicheweb.M), &(nicheweb.alpha), &(res.size));	

		if (checksum != 10 && checksum!=(int)(argc-1)/2) 	// Alles gesetzt?									
		 {	
			printf("Bitte gültige Eingabe für Parameter machen!\nProgramm wird beendet.\n");		
			return(0);		
		 }
	
//--Simulation---------------------------------------------------------------------------------------------------
	nicheweb.alpha = nicheweb.alpha/100;
	printf("alpha: %f\n", nicheweb.alpha);
	float Brech = (float)nicheweb.S/6;

	if(Brech>0) nicheweb.B = (int)(Brech + 0.5);
	
	else nicheweb.B =  (int)(Brech - 0.5);
	
	
	//int len	= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S);	// Länge des Rückabewerts

	gsl_vector *populationFIN 	= gsl_vector_calloc((nicheweb.Rnum + nicheweb.S)*(nicheweb.Y)*5 + (nicheweb.S) + 4*nicheweb.Y);				// Gleiche Länge wie Rückgabe von evolveNetwork
	gsl_vector *robustness		= gsl_vector_calloc(67);
	gsl_vector *robustnesstemp	= gsl_vector_calloc(67);
	gsl_vector *meanOfDataSqu	= gsl_vector_calloc(67);
	gsl_vector *meanSquOfData	= gsl_vector_calloc(67);
	gsl_vector *meanSquOfDatatemp	= gsl_vector_calloc(67);
	gsl_vector *standardDerivation	= gsl_vector_calloc(67);
	gsl_vector_set_zero(robustness);
	gsl_vector_set_zero(meanSquOfData);

	for(i = 0; i < L; i++)																							
	 { 			
		printf("\nStarte Durchlauf L = %i\n", i);
			
		nicheweb.network = SetNicheNetwork(nicheweb, res);
		populationFIN	 = EvolveNetwork(nicheweb);
		//printf("funDiv: %f\n",gsl_vector_get(populationFIN, 5*nicheweb.Y*(nicheweb.Rnum+nicheweb.S)+nicheweb.S+3*nicheweb.Y+3));
		// int j = 0;																					// Neues Netzwerk erzeugen
		// for(j=0; j<len; j++)printf("Netzwerk %i: %f", j, gsl_vector_get(nicheweb.network, j));
 	    // gsl_vector_add(populationFIN, EvolveNetwork(nicheweb));		
												
		gsl_vector_memcpy(robustnesstemp, EvaluateRobustness(populationFIN, nicheweb.Rnum, nicheweb.S, nicheweb.Y));	// Robustness Analyse

		gsl_vector_add(robustness,robustnesstemp);
		
		gsl_vector_memcpy(meanSquOfDatatemp,robustnesstemp);
		gsl_vector_mul(meanSquOfDatatemp,robustnesstemp);
		gsl_vector_add(meanSquOfData,meanSquOfDatatemp);
		
		printf("\nBeende Durchlauf L = %i\n", i);
	 }
	
	printf("L=%i\tspeciesini=%f\tspeciesfinal=%f\n", L, gsl_vector_get(robustness, 3)/L, gsl_vector_get(robustness, 9)/L);

//--Standardabweichung berechnen------------------------------------------------------------------------------------
	
	
	gsl_vector_memcpy(meanOfDataSqu,robustness);
	printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSqu,3));
	gsl_vector_mul(meanOfDataSqu,robustness);
	printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSqu,3));
	
	for(i =0; i<67; i++)
	{
	  gsl_vector_set(meanOfDataSqu, i, gsl_vector_get(meanOfDataSqu,i)/(L*L));
	  gsl_vector_set(meanSquOfData, i, gsl_vector_get(meanSquOfData,i)/L);
	  gsl_vector_set(standardDerivation, i, sqrt(gsl_vector_get(meanSquOfData,i)-gsl_vector_get(meanOfDataSqu,i)));
	}
	
	printf("S ist %f\n", gsl_vector_get(robustness,3));
	printf("Standardabweichung von S ist %f\n", gsl_vector_get(standardDerivation,3));
	printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSqu,3));
	printf("meanSquOfData ist %f\n", gsl_vector_get(meanSquOfData,3));
	
//--Abspeichern in File-------------------------------------------------------------------------------------	
	
	FILE *statistics;				// neuer FILE pointer
    char aims[255] = ORT;			// Ausgewähltes Verzeichnis
    char buffers[100];				// Speicher für Dateiname

    sprintf(buffers,"S%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%5.1falpha%5.2f.out",nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.T,L,res.size,nicheweb.alpha);		
	// sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

    statistics = fopen(strcat(aims, buffers),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
    // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

/*
    if(nicheweb.T==1&&nicheweb.Y==1&&nicheweb.d==0.0)

*/
      fprintf(statistics,"RSize\tS\tB\tM\tx\tY\tdpow\tT\tRob\tPerlok\tPerges\tSi_ges\tSi_TL1\tSi_TL2\tSi_TL3\tSi_TL4\tSi_TL>4\tSf_ges\tSf_TL1\tSf_TL2\tSf_TL3\tSf_TL4\tSf_TL>4\tBi_ges\tBi_TL1\tBi_TL2\tBi_TL3\tBi_TL4\tBi_TL>4\tBf_ges\tBf_TL1\tBf_TL2\tBf_TL3\tBf_TL4\tBf_TL>4\tSh_ges\tSh_TL1\tSh_TL2\tSh_TL3\tSh_TL4\tSh_TL>4\tBh_ges\tBh_TL1\tBh_TL2\tBh_TL3\tBh_TL4\tBh_TL>4\tSs_ges\tSs_TL1\tSs_TL2\tSs_TL3\tSs_TL4\tSs_TL>4\tBs_ges\tBs_TL1\tBs_TL2\tBs_TL3\tBs_TL4\tBs_TL>4\t1mit2\t2mit3\t3mit1\tFixp0\tFixp1\tFixp2\tFixp3\tFixp4\tFixp5\tFixp6\tFixp7\tRob2\tmetLoss\trPred\tiPred\tfDiv\talpha\n");

    fclose(statistics);
    statistics = fopen(aims,"a");												// fopen(*filename, a): schreibt am Ende der Datei weiter		

//--Daten in Datei schreiben---------------------------------------------------------------------------------------------------------------------------------				
    fprintf(statistics,"%5.1f\t%d\t%d\t%d\t%2.1f\t%d\t%2.1f\t%d\t", res.size, nicheweb.S, nicheweb.B, nicheweb.M, nicheweb.x, nicheweb.Y, nicheweb.d, nicheweb.T);		// Konsolenparameter im Namen

    for(i=0; i<51; i++)															// Rob, Perlok, Perges; +48 Elemente Spezies Informationen
      {
	gsl_vector_set(robustness, i, gsl_vector_get(robustness, i)/L);			// Scale (1/L)
        fprintf(statistics,"%5.3f\t", gsl_vector_get(robustness, i));			// %: Charakter; number.number: min Anzahl an Stellen die gedruckt werden; f: dezimal float
      }

	// muss hier skaliert werden oder nicht?
    for(i=51; i<62; i++)	//1mit2... Fixp1...7
      {
	fprintf(statistics,"%5.0f\t", gsl_vector_get(robustness, i));
      }

    for(i=62; i<67; i++)	//Rob2 -> regio 
      {
	gsl_vector_set(robustness, i, gsl_vector_get(robustness, i)/L);
        fprintf(statistics,"%5.3f\t", gsl_vector_get(robustness, i));
      }
    fprintf(statistics,"%5.3f\t", nicheweb.alpha);

    fprintf(statistics,"\n");
    
//--Standardabweichungen in die Datei (dritte Zeile) schreiben---------------------------------------------------------------    
    for(i= 0 ; i< 8; i++)
    {
      fprintf(statistics, "%d\t", 0);
    }
    
    for(i = 0 ; i<51; i++)
    {
      fprintf(statistics,"%5.1f\t", gsl_vector_get(standardDerivation, i));
    }
    
    for(i=51; i<62; i++)	//1mit2... Fixp1...7
    {
	fprintf(statistics,"n.b.\t");
    }

    for(i=62; i<67; i++)	//Rob2 -> regio 
    {
        fprintf(statistics,"%5.3f\t", gsl_vector_get(standardDerivation, i));
    }
    
    fprintf(statistics, "%d\t", 0);
    
    fprintf(statistics,"\n");
    
    fclose(statistics);															// Datei schließen

	printf("Simulation abgespeichert\n");
	
//--free----------------------------------------------------------------------------------------------------------------  
	free(nicheweb.network);
	
	gsl_vector_free(fixpunkte);
	gsl_vector_free(populationFIN);
	gsl_vector_free(robustness);	
	
	return(0);

}




















