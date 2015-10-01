/*
	
Diese Funktion wertet ein zeitlich entwickeltes Nischennetz in Bezug auf Robustheit aus. Der Rückgabewert enthält 61 Werte:
	Rob	Perlok	Perges
	Si_ges	Si_TL1	Si_TL2	Si_TL3	Si_TL4	Si_TL>4
	Sf_ges	Sf_TL1	Sf_TL2	Sf_TL3	Sf_TL4	Sf_TL>4
	Bi_ges	Bi_TL1	Bi_TL2	Bi_TL3	Bi_TL4	Bi_TL>4
	Bf_ges	Bf_TL1	Bf_TL2	Bf_TL3	Bf_TL4	Bf_TL>4
	Sh_ges	Sh_TL1	Sh_TL2	Sh_TL3	Sh_TL4	Sh_TL>4
	Bh_ges	Bh_TL1	Bh_TL2	Bh_TL3	Bh_TL4	Bh_TL>4
	Ss_ges	Ss_TL1	Ss_TL2	Ss_TL3	Ss_TL4	Ss_TL>4
	Bs_ges	Bs_TL1	Bs_TL2	Bs_TL3	Bs_TL4	Bs_TL>4
	1mit2	2mit3	3mit1
	Fixp0	Fixp1	Fixp2	Fixp3	Fixp4	Fixp5	Fixp6	Fixp7	Rob2
*/
#include "robustness.h"

#include <math.h>
#include <stdio.h>					
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector *EvaluateRobustness(gsl_vector* evolNetwork, int Rnum, int S, int Y)		//evolNetwork hat (Rnum+S)*Y*5+3+S Elemente
{
	gsl_vector *result = gsl_vector_calloc(68);

	printf("\nStarte Auswertung Robustness\n");

//--evolNetwork zerlegen-------------------------------------------------------------------------------------------

	gsl_vector_view yfin	= gsl_vector_subvector(evolNetwork, 0*Y*(Rnum+S), Y*(Rnum+S));					// Endergebnis
	gsl_vector_view y0	= gsl_vector_subvector(evolNetwork, 1*Y*(Rnum+S), Y*(Rnum+S));					// Startwerte der Populationsgrößen
	gsl_vector_view ymax	= gsl_vector_subvector(evolNetwork, 2*Y*(Rnum+S), Y*(Rnum+S));					// Maximalwerte nach t2
	gsl_vector_view ymin	= gsl_vector_subvector(evolNetwork, 3*Y*(Rnum+S), Y*(Rnum+S));					// Minimalwerte nach t2
	gsl_vector_view yavg	= gsl_vector_subvector(evolNetwork, 4*Y*(Rnum+S), Y*(Rnum+S));					// Durchschnitt nach t2
	gsl_vector_view tl	= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S), S);							// Troph. Level
  	gsl_vector_view metLoss	= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S)+S+3, Y);	
	gsl_vector_view resPred	= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S)+S+3+Y, Y);	
	gsl_vector_view intraPred= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S)+S+3+2*Y, Y);
  	gsl_vector_view funcDiversity= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S)+S+3+3*Y, Y);
	gsl_vector_view intraCompetition= gsl_vector_subvector(evolNetwork, 5*Y*(Rnum+S)+S+3+4*Y, Y);

	int i,l = 0;

	//for(l=0; l<Y; l++) printf("resPred=%f\n", gsl_vector_get(&resPred.vector, l));
	/*
	for(i=0; i<(Rnum+S)*Y; i++) printf("AVG=%f\n", gsl_vector_get(&yavg.vector, i));
	for(i=0; i< S; i++) printf("TL%i=%1.1f\n", i, gsl_vector_get(&tl.vector, i));
	*/
	

	int testf0 = (int)gsl_vector_get(evolNetwork, 5*(Rnum+S)+S+0);
	int testf1 = (int)gsl_vector_get(evolNetwork, 5*(Rnum+S)+S+1);
	int testf2 = (int)gsl_vector_get(evolNetwork, 5*(Rnum+S)+S+2);

	
	 
//--Fixpunktsuche im Ergebnis nach t2------------------------------------------------------------------------------
	
	printf("Fixpunktsuche\n");

	int fix3, fix4, fix5, fix6, fix7 = 0;
 
    for(i=0; i<(Rnum+S)*Y; i++)
    {
		float y_i	= gsl_vector_get(&yfin.vector, i);
		float max 	= gsl_vector_get(&ymax.vector, i);
		float min 	= gsl_vector_get(&ymin.vector, i);
				
      if(y_i <= 0)	//Spezies ausgestorben = Alle Werte erhöhen
      {
		fix3++;
		fix4++;
		fix5++;
		fix6++;
		fix7++;
      }
      else if ((max!=0) && (min!=0)) 	// zu keinem Zeitpunkt ausgestorben
      {
		if((((max-min)/(max+min))< 0.0001)	|| ((max-min)< 0.0001))	fix3++;		
		if(((max-min)/(max+min)) < 0.0001) 							fix4++;
		if((max-min)< 0.0001) 										fix5++;				// max Schwankung in Population kleiner als 0.0001
		if((((max-min)/(max+min))< 0.001)	|| ((max-min)< 0.001)) 	fix6++;
		if((((max-min)/(max+min))< 0.00001)	|| ((max-min)< 0.00001))fix7++;
      }
    }

//--Ergebnis Variablen---------------------------------------------------------------------------------------  
  double *speciesini	= (double *) calloc(6, sizeof(double));
  double *speciesfin	= (double *) calloc(6, sizeof(double));
  double *specieshubfin	= (double *) calloc(6, sizeof(double));
  double *speciessatfin	= (double *) calloc(6, sizeof(double));
  double *biomassini	= (double *) calloc(6, sizeof(double));
  double *biomassfin	= (double *) calloc(6, sizeof(double));
  double *biomasshubfin	= (double *) calloc(6, sizeof(double));
  double *biomasssatfin	= (double *) calloc(6, sizeof(double));
    
  double *regiorob		= (double *) calloc(S, sizeof(double));
  
  double Perlok	= 0;
  double Perges = 0;
  double Rob	= 0;
  double Rob2	= 0;

  int j, k	= 0;
  int rob2count	= 0;

  double metLossAvAll=0;
  double resPredAvAll=0;
  double intraPredAvAll=0;
  double funcDiversityAvAll=0;
  double intraCompetitionAvAll = 0;
  
//--Rob2 ausrechnen------------------------------------------------------------------------------------------------- 
  for(l=0; l<Y; l++)
    {
	  for(i= Rnum; i< Rnum+S; i++)
	 	{
	  	  regiorob[i-1] = regiorob[i-1] + gsl_vector_get(&yfin.vector, (Rnum+S)*l+i);	// Größe einer Spezies aus allen Patches addiert
		}
  	}
  
  for(i= 0; i<S; i++)
    {
	  if(regiorob[i] > 0) rob2count++;		// Spezies die überlebt haben (irgendwo)
    }
  
  Rob2 = (double)rob2count/(double)S;			// Anteil der Spezies die überlebt haben
    gsl_vector_set(result, 62, Rob2);			// Anteil der Spezies die auf irgendeinem Patch überlebt hat

	printf("Rob2=%f\n", Rob2);

//--Spezies und Biomassengrößen auswerten----------------------------------------------------------------------------
  for(l=0; l < Y; l++)
    {//*
		Rob		= 0;	// Für jeden Patch wieder nullen
		Perlok	= 0;

		for(i=0;i<6;i++)
			{
			  speciesfin[i]=0;
			  biomassfin[i]=0;
			  specieshubfin[i]=0;
			  biomasshubfin[i]=0;
			  speciessatfin[i]=0;
			  biomasssatfin[i]=0;
			  speciesini[i]=0;
			  biomassini[i]=0;
			}

		//printf("species ini = %f\n", speciesini[0]);
	
		for(i = Rnum; i< (Rnum+S); i++)			 	// i die Spezies durch
		  {
			  float y	= gsl_vector_get(&yfin.vector, l*(Rnum+S)+i);		// jeweils die Vektorelemente
			  float avg	= gsl_vector_get(&yavg.vector, l*(Rnum+S)+i);			
			  float y0i = gsl_vector_get(&y0.vector,   l*(Rnum+S)+i);
			//  int tempi	= (int)fmod(i, (Rnum+S));						// tempi zählt die Spezies durch (fmod gibt Rest zurück Bsp: 2:3 -> Rest 2, da 0*3+2 = 2)
			  
			  speciesini[0]	= speciesini[0] + 1;							// Anfangsspezies Anzahl == S
			  biomassini[0]	= biomassini[0] + y0i;							// Biomasse am Anfang 	 == y0 
			  
			  k = (int)gsl_vector_get(&tl.vector, i-Rnum);					// tropl
				//printf("k=%i\n", k);
			  if (k>4)
				k = 5;
				
				//printf("k=%i\n", k);

			  if(k!=0){

			  speciesini[k]	= speciesini[k] + 1;			// Spezies nach TL sortiert (wie viele pro Level)
			  biomassini[k]	= biomassini[k] + y0i;			// Biomasse pro TL

				}
			//--Falls Spezies überlebt----------------------------------------------------------------------------------
			  if(y !=0)																
			  {
				speciesfin[0] = speciesfin[0] + 1;			// Anzahl Spezies, die überlebt
				biomassfin[0] = biomassfin[0] + avg;		// Durchschnittliche Biomasse der Spezies wird mit berücksichtigt
					
				if(k!=0){

				speciesfin[k] = speciesfin[k] + 1;									
				biomassfin[k] = biomassfin[k] + avg;

				}

				if(l==0)	// Patch1
					{
					  biomasshubfin[0] = biomasshubfin[0] + avg;	// biomasshubfin(0) = Summe(avg_i) Y=0
					  specieshubfin[0] = specieshubfin[0] + 1;		// specieshubfin(0) = S(überlebt)
					  
					  biomasshubfin[k] = biomasshubfin[k] + avg;	// biomasshubfin(k) = Summe(avg_tl)
					  specieshubfin[k] = specieshubfin[k] + 1;		// specieshubfin(k) = Spezies mt TL=k (überlebt)
					}
	
				if(l==1)	// Patch2
					{
					  biomasssatfin[0] = biomasssatfin[0] + avg;	// biomasssatfin(0) = Summe(avg_i) Y=1
					  speciessatfin[0] = speciessatfin[0] + 1;		// biomasssatfin(0) = S(überlebt)
					  
					  biomasssatfin[k] = biomasssatfin[k] + avg;	// biomasssatfin(k) = Summe(avg_tl) Y=1 
					  speciessatfin[k] = speciessatfin[k] + 1; 		// biomasssatfin(k) = S mit TL = k(überlebt)
					}
			  }
			//------------------------------------------------------------------------------------------------------------
			
			//printf("species fin = %f\n", speciesfin[0]);

		}//Ende Spezies Anzahl Analyse

			
		//printf("Überlebende Spezies ausgewertet\n");

		
		

		metLossAvAll += gsl_vector_get(&metLoss.vector, l);
		resPredAvAll += gsl_vector_get(&resPred.vector, l);
		intraPredAvAll += gsl_vector_get(&intraPred.vector, l);
		funcDiversityAvAll += gsl_vector_get(&funcDiversity.vector, l);
		intraCompetitionAvAll += gsl_vector_get(&intraCompetition.vector,l);
		  //printf("funcDiversity schleife: %f\n",funcDiversityAvAll);
		//printf("metLossAvAll: %f\n",metLossAvAll);
		//printf("resPredAvAll: %f\n",resPredAvAll);
		//--speichern in result pro schleifenaufruf---------------------------------------------------------------

		Rob 	= speciesfin[0]/speciesini[0];												// Wenn alle Spezies überleben ist Rob == 1 (für den aktuellen Patch!)
		printf("Patch %i species ini = %f\tspecies fin = %f\n", l, speciesini[0], speciesfin[0]);
		printf("Rob = %f\n", Rob);
		gsl_vector_set(result, 0, gsl_vector_get(result, 0) + Rob);      	//result(0) = Rob Addition über alle Patches
		
		if(Rob == 1.0) Perlok=1.0;							// Da pro Patch ausgewertet wird ist Perlok = Rob (lokale Persistence)		
		gsl_vector_set(result, 1, gsl_vector_get(result, 1) + Perlok);
		    
		Perges	= Perges + speciesfin[0];					// Summe Spezies die überlebt alle Patches
		
		//	printf("species ini = %f\n", speciesini[0]);		  

		for(i=0; i<6; i++)
		  {
			gsl_vector_set(result, 3 +i, gsl_vector_get(result, 3 + i) + speciesini[i]);
			gsl_vector_set(result, 9 +i, gsl_vector_get(result, 9 + i) + speciesfin[i]);
			gsl_vector_set(result, 15+i, gsl_vector_get(result, 15+ i) + biomassini[i]);
			gsl_vector_set(result, 21+i, gsl_vector_get(result, 21+ i) + biomassfin[i]);
			gsl_vector_set(result, 27+i, gsl_vector_get(result, 27+ i) + specieshubfin[i]);
			gsl_vector_set(result, 33+i, gsl_vector_get(result, 33+ i) + biomasshubfin[i]);
			gsl_vector_set(result, 39+i, gsl_vector_get(result, 39+ i) + speciessatfin[i]);
			gsl_vector_set(result, 45+i, gsl_vector_get(result, 45+ i) + biomasssatfin[i]);
		  }
		
    }//*


    
    
//--result skalieren ---------------------------------------------------------------------------------------------------
 gsl_vector_set(result, 0, gsl_vector_get(result, 0)/(double)Y);		// Rob skalieren auf pro Patch
 gsl_vector_set(result, 1, gsl_vector_get(result, 1)/(double)Y);		// Perges ebenso

  if(Perges == (double)S*(double)Y)										// Falls alles überlebt hat is result(2) = 1
	gsl_vector_set(result, 2, 1.0);

  for(i=3; i<27; i++)													// Skaliere speciesini - biomassfin auf 1 Patch
	gsl_vector_set(result, i, gsl_vector_get(result, i)/(double)Y);

//--Fixvariablen speichern-------------------------------------------------------------------------------------  
  if(testf0==1)	gsl_vector_set(result, 54, 1);
     	   else gsl_vector_set(result, 54, 0);

  if(testf1==1)	gsl_vector_set(result, 55, 1);
     	   else gsl_vector_set(result, 55, 0);

  if(testf2==1)	gsl_vector_set(result, 56, 1);
     	   else gsl_vector_set(result, 56, 0);

  if(fix3==(Rnum+S)*Y)	gsl_vector_set(result, 57, 1);
     	   		   else gsl_vector_set(result, 57, 0);

  if(fix4==(Rnum+S)*Y)	gsl_vector_set(result, 58, 1);
     	   		   else gsl_vector_set(result, 58, 0);

  if(fix5==(Rnum+S)*Y)	gsl_vector_set(result, 59, 1);
     	   		   else gsl_vector_set(result, 59, 0);

  if(fix6==(Rnum+S)*Y)	gsl_vector_set(result, 60, 1);
     	   		   else gsl_vector_set(result, 60, 0);

  if(fix7==(Rnum+S)*Y)	gsl_vector_set(result, 61, 1);
     	   		   else gsl_vector_set(result, 61, 0);

//--Result 51, 52, 53- Synchronisation?-----------------------------------------------------------------------------------------------------
  l=0;
  
  if(Y==1)	// Für 1-Patch Fall Was bringt diese Analyse?!?
  {
    for(i=0; i<Rnum+S; i++)
      {
		float y_i = gsl_vector_get(&yfin.vector, i);

        if((y_i <= y_i + 0.00001) && (y_i >= y_i -0.00001)) 	// Falls y auf 0.00001 genau ist zähle l hoch
		l++;
      }

    if(l == (Rnum+S)) gsl_vector_set(result, 51, 1);			// Alle y im Intervall
    else gsl_vector_set(result, 51, 0);							// mind 1 y Wert nicht im Intervall
  }

  else if(Y==2)	// Für 2-Patch Fall
  {
    for(i=0; i<Rnum+S; i++)
      {
	    float y_i = gsl_vector_get(&yfin.vector, i);			// Population in Patch 1
		float y_l = gsl_vector_get(&yfin.vector, (Rnum+S)+i);	// Population in Patch 2

        if((y_i <= y_l +0.00001) && (y_i >= y_l - 0.00001)) 	// Falls Wert von Spezies S in Patch 1 fast genau Wert von S in Patch 2 ist 
		l++;													// Erhöhe l
      }

	if(l == (Rnum+S)) gsl_vector_set(result, 51, 1);			// Alle Spezies in 1 und 2 fast gleiche Population
    else gsl_vector_set(result, 51, 0);							// mind 1 Spezies nicht die gleiche Population
  }

  else			// mehr als 2 Patches
  {
    for(i=0; i<Rnum+S; i++)
      {
		float y_i = gsl_vector_get(&yfin.vector, 	0 *(Rnum+S)+i);			// Population in Patch 1
		float y_l = gsl_vector_get(&yfin.vector, 	1 *(Rnum+S)+i);			// Population in Patch 2
		float y_j = gsl_vector_get(&yfin.vector, 	2 *(Rnum+S)+i);			// Population in Patch 3
		float y_k = gsl_vector_get(&yfin.vector, (Y-1)*(Rnum+S)+i);		// Population im letzten Patch
		
	    if((y_i <= y_l + 0.00001) && (y_i >= y_l - 0.00001)) 			// Population in 1 fast wie in 2	
			l++;
        if((y_j <= y_l + 0.00001) && (y_j >= y_l - 0.00001)) 			// Population in 2 fast wie in 3
			j++;
	    if((y_j <= y_k + 0.00001) && (y_j >= y_k - 0.00001))			// Population in 3 fast wie im letzten Patch
			k++;
      }

	 if(l == (Rnum+S)) gsl_vector_set(result, 51, 1);			// 1 mit 2 
		else gsl_vector_set(result, 51, 0);							
	 if(j == (Rnum+S)) gsl_vector_set(result, 52, 1);			// 2 mit 3
		else gsl_vector_set(result, 52, 0);							
	 if(k == (Rnum+S)) gsl_vector_set(result, 53, 1);			// 3 mit letztem
		else gsl_vector_set(result, 53, 0);							
  }
  gsl_vector_set(result, 63, metLossAvAll/Y);
  gsl_vector_set(result, 64, resPredAvAll/Y);
  gsl_vector_set(result, 65, intraPredAvAll/Y);
  gsl_vector_set(result, 66, funcDiversityAvAll/Y);
  gsl_vector_set(result, 67, intraCompetitionAvAll/Y);

//--free----------------------------------------------------------------------------------------------------------------------------  

  free(speciesfin);
  free(speciesini);
  free(biomassfin);
  free(biomassini);
  free(specieshubfin);
  free(speciessatfin);
  free(biomasshubfin);
  free(biomasssatfin);
  free(regiorob);
  
//--return---------------------------------------------------------------------------------------------------------------------------

  return result;

}//Ende Robustness Analyse

