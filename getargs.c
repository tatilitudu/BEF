/*	Definition von getArgs
 	Liest die Parameter der Simulation aus der Kommandozeile und wandelt sie in numerische Werte um
 	Die Eingabe erfolgt über -Param_Leerzeichen_WertVonParam_Leerzeichen (ohne Unterstriche),
	wobei die Leerzeichen als Zeilentrenner dienen.

  	-> in Felder Darstellung dann für die Beispiel-Eingabe
	
	 -S 1 -B 1 -T 0 -d -7 -L 1 -Y 1 -x 0.0 -M "w"

	[1][0] - [1][1] S
	[2][0] 1
	[3][0] - [3][1] B
	[4][0] 1  
	[5][0] - [5][1] T
	[6][0] 0
	[7][0] - [7][1] d
	[8][0] - [8][1] 7
	[9][0] - [9][1] L
	[10][0] 1  
	[11][0] - [11][1] Y
	[12][0] 1
	[13][0] - [13][1] x
	[14][0] 0 [14][1] . [14][2] 0
	[15][0] - [15][1] M  
	[16][0] k

	Im Gegensatz zur Version von Sebastian get_args.c werden hier keine globalen Variablen verwendet. 
	Die Funktion getArgs() gibt ein Array zurück, in denen die Konsolenparameter stehen. Wenn keine Eingabe gemacht wurde werden
	die default Parameter gesetzt und eine entsprechende Warnung ausgegeben. KEINE AHNUNG OB DAS WAS NUTZT
  
*/
#include <string.h>						// string modification functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard

int getArgs(int argc, char** argv, int* S_value, int* B_value, int* T_value, double* d_value, int* L_value, int* Y_value, double* x_value, int* M_value, double* alpha_value, double* R_value)
{
  int i = 0;

  // Initialisiere Check-Variablen
  int S_check	= 0;
  int B_check	= 0;
  int T_check	= 0;
  int d_check	= 0;
  int L_check	= 0;
  int Y_check	= 0;
  int x_check	= 0;
  int M_check	= 0;
  int R_check 	= 0;
  int alpha_check = 0;

  int checksum	= 0;		// Wird am Ende ausgegeben

  // Start bei [1][i], da in argv[0][i] der Name des Programmes steht 
  // Solange i kleiner als der Argument Counter (Länge = Anzahl Parameter+2) sind noch Paramter einzulesen
  for (i = 1; i < argc; i++) 	
  {
	
  /* argv[i][j] enthält die über die Konsole eingegebenen Parameter, siehe Erklärung oben.
     atoi(*char) bzw. atof(*char) (ASCII to Int/Float) machen aus einer gegebenen Zeichenkette dessen numerischen Wert (int oder float).
  */

    if (argv[i][0] == '-') 	// "-" als Indikator für neuen Parameter	
    {
      switch (argv[i][1]) 	// Ab Feld [i][1] können die Parameter identifiziert werden, der Wert befindet sich ein Feld weiter in argv[++i]
      {
		case 'S':	*S_value = atoi(argv[++i]); // Der Speicher auf den S_value zeigt wird mit "atoi(argv[++i])" = dem numerischen Werte des nächsten Feldes, beschrieben
					S_check  = 1;
					printf("S gesetzt auf %i \n", *S_value);
					break;
			
		case 'B':	*B_value = atoi(argv[++i]); // Der Speicher auf den S_value zeigt wird mit "atoi(argv[++i])" = dem numerischen Werte des nächsten Feldes, beschrieben
					B_check  = 1;
					printf("B gesetzt auf %i \n", *B_value);
					break;

		case 'T':	*T_value = atoi(argv[++i]);
					T_check  = 1;
					printf("T gesetzt auf %i \n", *T_value);
					break;
	
		case 'd':	*d_value = atof(argv[++i]);
					d_check  = 1;
					printf("d gesetzt auf %f \n", *d_value);
					break;
			
		case 'L':	*L_value = atoi(argv[++i]);
					L_check  = 1;
					printf("L gesetzt auf %i \n", *L_value);
					break;
			
		case 'Y':	*Y_value = atoi(argv[++i]);
					Y_check  = 1;
					printf("Y gesetzt auf %i \n", *Y_value);
					break;
	
		case 'x':	*x_value = atof(argv[++i]);
					x_check  = 1;
					printf("x gesetzt auf %f \n", *x_value);
					break;
	
		case 'M':	*M_value = atoi(argv[++i]);
					M_check  = 1;
					printf("M gesetzt auf %i \n", *M_value);
					break;
					
		case 'a':	*alpha_value = atof(argv[++i]);
					alpha_check = 1;
					printf("alpha vorübergehend gesetzt auf %f \n", *alpha_value);
					break;

		case 'R':	*R_value = atof(argv[++i]);
					R_check  = 1;
					printf("Resource Größe ist %f \n", *R_value);
					break;
	


		default:	fprintf(stderr,"Unknown switch: %s\n", argv[i]); // unbekannter Parametername -> Fehler in Konsole
      }
    }
  }

    if (S_check == 0) printf("Missing switch: -S\n");
    if (B_check == 0) printf("Missing switch: -B\n");
    if (T_check == 0) printf("Missing switch: -T\n");
    if (d_check == 0) printf("Missing switch: -d\n");
    if (L_check == 0) printf("Missing switch: -L\n");
    if (Y_check == 0) printf("Missing switch: -Y\n");
    if (x_check == 0) printf("Missing switch: -x\n");
    if (M_check == 0) printf("Missing switch: -M\n");
    if (R_check == 0) printf("Missing switch: -R\n");
    if (alpha_check == 0) printf("Missing switch: -a\n");
    
    checksum = S_check + B_check + T_check + d_check + L_check + Y_check + x_check + M_check + R_check + alpha_check;
	
	printf("checksum ist %i \n", checksum);
    
    return checksum;	// Gibt Prüfsumme zurück, diese sollte 8 sein
}

