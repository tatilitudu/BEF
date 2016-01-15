/*
Dieses Programm erzeugt eine Matrix, die Y verbundene Lebensräume (Patches) simuliert. 
Ein Eintrag 1 bedeutet, man kann von Patch "Zeile" nach Patch "Spalte" gelangen. 0 bedeutet es gibt zwischen Zeile und Spalte keine Verbindung. 
Mögliche Verbindungsmuster sind 

0: Keine Verbindung -> Matrix enthält nur Null als Einträge
1: Kette -> Nebendiagonalen sind Einsen
2: Ring -> Nebendiagonalen sind Einsen, sowie die Einträge [0][Y] und [Y][0]
3: Vollständig verbunden -> Matrix enthält nur Einsen als Einträge, bis auf der Diagonalen
4: Stern -> Erste Zeile und erste Spalte enthalten Einsen, bis auf [0][0]

*/

#include "topology.h"				// deklariert im Topology Header
#include <gsl/gsl_blas.h>
#include <math.h>					// math functions
#include <stdio.h>					// output functions
#include <stdlib.h>					// standard
#include <gsl/gsl_vector.h>

gsl_matrix *SetTopology(int Y, int T, gsl_matrix* D){

int i = 0;



	switch(T)
	{
		//--Keine-----------------------------------------------------
			case 0:	 gsl_matrix_set_zero(D);
					 printf("Migration ausgeschaltet\n");	
					 break;	

		//--Kette------------------------------------------------------
			case 1:  for(i=0; i<Y; i++)						
					  {
						if(i!=Y-1)	gsl_matrix_set(D,i,i+1,1.);
					    if(i!=0)	gsl_matrix_set(D,i,i-1,1.);
					  }
					 gsl_matrix_set(D, 0, 0, 0.);
					 printf("Topologie auf Kette gesetzt\n");
					 break;	

		//--Ring-------------------------------------------------------
			case 2:  for(i=0; i<Y; i++)									
					  {
						 if(i!= Y-1) gsl_matrix_set(D, i, i+1, 1.);
						 else 		 gsl_matrix_set(D, i, 0, 1.);

						 if(i!= 0)	 gsl_matrix_set(D, i, i-1, 1.);
						 else		 gsl_matrix_set(D, i, Y-1, 1.);
					  }
					 gsl_matrix_set(D, 0, 0, 0.);
					 printf("Topologie auf Ring gesetzt\n");
					 break;	

		//--Total------------------------------------------------------
			case 3:  gsl_matrix_set_all(D,1.);
   					 for(i=0; i<Y; i++) 
					 gsl_matrix_set(D, i, i, 0.);
					 printf("Topologie auf total gesetzt\n");
					 break;	

		//--Stern------------------------------------------------------
			case 4:	 for(i=0; i<Y; i++) 
					 {
						gsl_matrix_set(D, 0, i, 1.);
						gsl_matrix_set(D, i, 0, 1.);
					 }
					 gsl_matrix_set(D, 0, 0, 0.);
					 printf("Topologie auf Stern gesetzt\n");
					 break;	

		//--Gitter------------------------------------------------------
			/*case 5:	 for(i=0; i<Y; i++) 
					 {
						if(i!==0 || i!= Y || mod(i, sqrt(Y)-1)!=0) 
								
							{	gsl_matrix_set(D, i,   i+1, 1.);
								gsl_matrix_set(D, i,   i+1, 1.);
								gsl_matrix_set(D, i+1, i,   1.);
								gsl_matrix_set(D, i-1, i,   1.);
							}	
						else if(i==0) 
							{	gsl_matrix_set(D, i,   i+1, 1.);
								gsl_matrix_set(D, i,   i+1, 1.);
								gsl_matrix_set(D, i+1, i,   1.);
								gsl_matrix_set(D, i-1, i,   1.);
							}	


gsl_matrix_set(D, 0, i, 1.);			
						
						gsl_matrix_set(D,i,0,1.);
					 }
					 gsl_matrix_set(D,0,0,0.);
					 printf("Topologie auf Gitter gesetzt\n");
					 break;	

*/

		//--Default-----------------------------------------------------
			default: printf("Unbekannte Eingabe für T, Migrationsmatrix leer!");
					 break;
				
	}

  return D;

}

