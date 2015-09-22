
#include "structs.h"

#include <math.h>					// math functions
#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>


double* predOnRes(struct foodweb nicheweb, const double y[], double* resPred) 
{

  int S 	 	= nicheweb.S;
  int Y 	     	= nicheweb.Y;
  int Rnum 		= nicheweb.Rnum;
  gsl_vector *network 	= nicheweb.network;						// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
  double lambda		= nicheweb.lambda;
  double aij		= nicheweb.aij;
  double hand		= nicheweb.hand;
  
  int i,j,l;

  /* Massen rausholen */
  gsl_vector_view A_view = gsl_vector_subvector(network, 0, (Rnum+S)*(Rnum+S));						// Fressmatrix A als Vektor
  gsl_matrix_view EA_mat = gsl_matrix_view_vector(&A_view.vector, (Rnum+S), (Rnum+S));				// A als Matrix_view
  gsl_matrix *EAmat	   = &EA_mat.matrix;		// A als Matrix

  gsl_vector_view M_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
  gsl_vector *Mvec	   = &M_vec.vector;

  
  double ytemp[(Rnum+S)*Y];		// tempvector for populations and efforts
  for(i=0;i<(Rnum+S)*Y;i++)
    ytemp[i]=y[i];





  /* Auslesen von ytemp = y[]; sind Population */
  gsl_vector_view yfd_vec=gsl_vector_view_array(ytemp,(Rnum+S)*Y);
  gsl_vector *yfdvec=&yfd_vec.vector;				// populations and efforts for later use
  
  gsl_matrix *EAmat_new= gsl_matrix_calloc(Rnum+S, Rnum+S);
  
  gsl_matrix_memcpy(EAmat_new,EAmat);
  
  for(i=0;i<Rnum+S;i++)
  {
    for(j=Rnum;j<Rnum+S;j++)
    {
      gsl_matrix_set(EAmat_new,i,j,0);
    }
  }
  

  for(i=Rnum;i<Rnum+S;i++)
  {
    for(j=0;j<Rnum+S;j++)
    {
      if(gsl_matrix_get(EAmat_new,i,j)==1)
      {
	//printf("%3.1f\t", gsl_vector_get(network, (S+Rnum)*(S+Rnum)+1+(Y*Y)+1+(Rnum+S)+(j-Rnum)));
      }
    }
  }
  
  /* Initialisierungen */
  gsl_matrix *AFgsl=gsl_matrix_calloc(Rnum+S, Rnum+S);		// matrix of foraging efforts
  
  gsl_matrix *Emat=gsl_matrix_calloc(Rnum+S, Rnum+S);		// gsl objects for calculations of populations 
  gsl_vector *tvec=gsl_vector_calloc(Rnum+S);
  gsl_vector *rvec=gsl_vector_calloc(Rnum+S);
  gsl_vector *svec=gsl_vector_calloc(Rnum+S);
  gsl_vector *resPredTemp=gsl_vector_calloc(Rnum+S);
  
  for(l=0;l<Y;l++)						// start of patch solving
  {
    /* Initialisierungen */
    gsl_matrix_set_zero(AFgsl);					// reset gsl objects for every patch
    gsl_matrix_set_zero(Emat);
    gsl_vector_set_zero(tvec);
    gsl_vector_set_zero(rvec);
    gsl_vector_set_zero(svec);
    
    
    /* Je Vektoren von (Res+S) Elementen */

    /* yfdvec enthält die Population */
    gsl_vector_view y_vec=gsl_vector_subvector(yfdvec,(Rnum+S)*l,(Rnum+S));
    gsl_vector *yvecpre=&y_vec.vector;
    
    /* Kopie von EAmat erstellen */
    gsl_matrix_memcpy(AFgsl,EAmat_new);

    for(i=0;i<Rnum+S;i++)
    {
      /* Nehme i-te Zeile aus A */
      gsl_vector_view tempp=gsl_matrix_row(AFgsl,i);
      /* Summiere Absolutwerte der Zeile */
      double temp1;	
      temp1=gsl_blas_dasum(&tempp.vector);
      if(temp1!=0)
      {
	/* Teile die Einträge, die nicht- Null sind durch Anzahl an nicht-Nullen in dieser Zeile*/ 
	/* und setzte diesen Wert dann an den entsprechenden Platz */
	/* Man erhält also eine prozentuale Verbindung */
	for(j=0;j<Rnum+S;j++)
	  gsl_matrix_set(AFgsl,i,j,(gsl_matrix_get(AFgsl,i,j)/temp1));
      }
    }
  
  /* aij ist Attackrate; AFgsl ist jetzt normiert- also fij  */
    gsl_matrix_memcpy(Emat,EAmat_new);
    gsl_matrix_scale(Emat,aij);					//  Emat(i,j) = a(i,j)
    gsl_matrix_mul_elements(Emat,AFgsl);			//  Emat(i,j) = a(i,j)*f(i,j)

    
    /*  hand =  handling time */
    /* Berechnung wie aus Paper */
    gsl_vector_memcpy(svec,yvecpre);				// s(i)=y(i)
    gsl_vector_scale(svec, hand);				// s(i)=y(i)*h
    gsl_blas_dgemv(CblasNoTrans,1,Emat,svec,0,rvec);		// r(i)=Sum_k h*a(i,k)*f(i,k)*y(k)
    gsl_vector_add_constant(rvec,1);				// r(i)=1+Sum_k h*a(i,k)*f(i,k)*y(k)
    
    gsl_vector_memcpy(tvec,Mvec);				// t(i)=masse(i)^(-0.25)
    gsl_vector_div(tvec,rvec);					// t(i)=masse(i)^(-0.25)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))
    gsl_vector_mul(tvec,yvecpre);					// t(i)=masse(i)^(-0.25)*y(i)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))

  
  //gsl_vector_memcpy(ydottemp,ydotvec);

//   for(i=Rnum; i<Rnum+S; i++)
//   {
//     /* Benötige nur Resource */		
//     gsl_vector_set(yvecpre,i,0);
//   }
  //printf("yvecpre: %f\n",gsl_vector_get(yvecpre,0));
  
  gsl_blas_dgemv(CblasNoTrans,lambda,Emat,yvecpre,0,resPredTemp);	// ydot(i) = Sum_j lambda*a(i,j)*f(i,j)*y(j)
  gsl_vector_mul(resPredTemp,tvec);				// ydot(i) = Sum_j lambda*a(i,j)*f(i,j)*y(j)*t(i)
  resPred[l] = gsl_blas_dasum(resPredTemp);			// PredOnRes = Sum_i ydot(i)
   //printf("resPred %f\n",resPred[0]);   
  }
  
 gsl_matrix_free(Emat);  
 gsl_matrix_free(EAmat_new);
 gsl_matrix_free(AFgsl); 
 gsl_vector_free(tvec);
 gsl_vector_free(rvec);
 gsl_vector_free(svec);
 gsl_vector_free(resPredTemp);
 
 return 0;
 
}


double* intraguildPred(struct foodweb nicheweb, const double y[], double* intraPred)
{
  int i,j,l;

  int S 	 	= nicheweb.S;
  int Y 	     	= nicheweb.Y;
  int Rnum 		= nicheweb.Rnum;
  gsl_vector *network 	= nicheweb.network;						// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
  double lambda		= nicheweb.lambda;
  double aij		= nicheweb.aij;
  double hand		= nicheweb.hand;

  /* Massen rausholen */
  gsl_vector_view A_view = gsl_vector_subvector(network, 0, (Rnum+S)*(Rnum+S));						// Fressmatrix A als Vektor
  gsl_matrix_view EA_mat = gsl_matrix_view_vector(&A_view.vector, (Rnum+S), (Rnum+S));				// A als Matrix_view
  gsl_matrix *EAmat	   = &EA_mat.matrix;		// A als Matrix

  gsl_vector_view M_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
  gsl_vector *Mvec	   = &M_vec.vector;				// massvector: M(i)=m^(-0.25)
  
  double ytemp[(Rnum+S)*Y];		// tempvector for populations and efforts
  for(i=0;i<(Rnum+S)*Y;i++)
    ytemp[i]=y[i];

  /* Alles view_array */
  
  /* Auslesen von ytemp = y[]; sind Population */
  gsl_vector_view yfd_vec=gsl_vector_view_array(ytemp,(Rnum+S)*Y);
  gsl_vector *yfdvec=&yfd_vec.vector;				// populations and efforts for later use
  
 
  
  
  /* Initialisierungen */
  gsl_matrix *AFgsl=gsl_matrix_calloc(Rnum+S, Rnum+S);		// matrix of foraging efforts
  
  gsl_matrix *Emat=gsl_matrix_calloc(Rnum+S, Rnum+S);		// gsl objects for calculations of populations 
  gsl_vector *tvec=gsl_vector_calloc(Rnum+S);
  gsl_vector *rvec=gsl_vector_calloc(Rnum+S);
  gsl_vector *svec=gsl_vector_calloc(Rnum+S);
  gsl_vector *yvecint=gsl_vector_calloc(Rnum+S);
  gsl_vector *intraPredTemp=gsl_vector_calloc(Rnum+S);
  
  
  for(l=0;l<Y;l++)						// start of patch solving
  {
    /* Initialisierungen */
    gsl_matrix_set_zero(AFgsl);					// reset gsl objects for every patch
    gsl_matrix_set_zero(Emat);
    gsl_vector_set_zero(tvec);
    gsl_vector_set_zero(rvec);
    gsl_vector_set_zero(svec);
   
    
    /* Je Vektoren von (Res+S) Elementen */


    /* yfdvec enthält die Population */
    gsl_vector_view y_vec=gsl_vector_subvector(yfdvec,(Rnum+S)*l,(Rnum+S));
    gsl_vector *yveci=&y_vec.vector;
    
    gsl_vector_memcpy(yvecint,yveci);
    
    /* Kopie von EAmat erstellen */
    gsl_matrix_memcpy(AFgsl,EAmat);

    for(i=0;i<Rnum+S;i++)
    {
      /* Nehme i-te Zeile aus A */
      gsl_vector_view tempp=gsl_matrix_row(AFgsl,i);
      
      /* Summiere Absolutwerte der Zeile */
      double temp1;	
      temp1=gsl_blas_dasum(&tempp.vector);
      if(temp1!=0)
      {
	/* Teile die Einträge, die nicht- Null sind durch Anzahl an nicht-Nullen in dieser Zeile*/ 
	/* und setzte diesen Wert dann an den entsprechenden Platz */
	/* Man erhält also eine prozentuale Verbindung */
	for(j=0;j<Rnum+S;j++)
	  gsl_matrix_set(AFgsl,i,j,(gsl_matrix_get(AFgsl,i,j)/temp1));
      }
    }
  
  /* aij ist Attackrate; AFgsl ist jetzt normiert- also fij  */
    gsl_matrix_memcpy(Emat,EAmat);
    gsl_matrix_scale(Emat,aij);					//  Emat(i,j) = a(i,j)
    gsl_matrix_mul_elements(Emat,AFgsl);			//  Emat(i,j) = a(i,j)*f(i,j)

    
    /*  hand =  handling time */
    /* Berechnung wie aus Paper */
    gsl_vector_set(yvecint,0,0);
    //printf("y: %f\n",gsl_vector_get(yvecint,0));
    gsl_vector_memcpy(svec,yvecint);				// s(i)=y(i)
    gsl_vector_scale(svec, hand);				// s(i)=y(i)*h
    gsl_blas_dgemv(CblasNoTrans,1,Emat,svec,0,rvec);		// r(i)=Sum_k h*a(i,k)*f(i,k)*y(k)
    gsl_vector_add_constant(rvec,1);				// r(i)=1+Sum_k h*a(i,k)*f(i,k)*y(k)
    
    gsl_vector_memcpy(tvec,Mvec);				// t(i)=masse(i)^(-0.25)
    gsl_vector_div(tvec,rvec);					// t(i)=masse(i)^(-0.25)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))
    gsl_vector_mul(tvec,yvecint);				// t(i)=masse(i)^(-0.25)*y(i)/(1+Sum_k h*a(i,k)*f(i,k)*y(k))

    gsl_blas_dgemv(CblasNoTrans,lambda,Emat,yvecint,0,intraPredTemp);	// intraPredTemp(i)=Sum_j lambda*a(i,j)*f(i,j)*y(j)
    gsl_vector_mul(intraPredTemp,tvec);
    
    
//    double asum =0;
//     for(i=0;i<S+Rnum;i++)
//     {
//       for(j=0;j<Rnum+S;j++)
//       {
// 	asum += gsl_matrix_get(EAmat,i,j);
// 	//printf("Element i: %i, j: %i von A: %f\n",i,j,gsl_matrix_get(EAmat,i,j));
//       }
//       //printf("intraPredTemp %f \t i %i\n",gsl_vector_get(intraPredTemp,i),i);
//     }
//     printf("Summe aller Elemente von A: %f\n",asum);
    
    intraPred[l] = gsl_blas_dasum(intraPredTemp);
    for(i=0;i<S+Rnum;i++)
    {
      //printf("y %f \t i %i\n",gsl_vector_get(yvecint,i),i);
    }
    //printf("intraguildPred %f\n",intraPred[0]);
  }
  /* Speicher befreien */
  gsl_matrix_free(Emat); 
  gsl_matrix_free(AFgsl);  
  
  gsl_vector_free(tvec);
  gsl_vector_free(rvec);
  gsl_vector_free(svec);
  gsl_vector_free(intraPredTemp);
  
  return 0;
}


double* metabolicLoss(struct foodweb nicheweb, const double y[], double* metLoss)
{

  int S 	 	= nicheweb.S;
  int Y 	     	= nicheweb.Y;
  int Rnum 		= nicheweb.Rnum;
  double alpha		= nicheweb.alpha;
  gsl_vector *network 	= nicheweb.network;						// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
  
  int i,l;


  /* Massen rausholen */
  gsl_vector_view M_vec  = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
  gsl_vector *Mvec	   = &M_vec.vector;
  
  double ytemp[(Rnum+S)*Y];		// tempvector for populations and efforts
  for(i=0;i<(Rnum+S)*Y;i++)
    ytemp[i]=y[i];

  /* Alles view_array */
  
 
  /* Auslesen von ytemp = y[]; sind Population */
  gsl_vector_view yfd_vec=gsl_vector_view_array(ytemp,(Rnum+S)*Y);
  gsl_vector *yfdvec=&yfd_vec.vector;				// populations and efforts for later use
  
  
  
  /* Initialisierungen */
  gsl_vector *svec=gsl_vector_calloc(Rnum+S);
  gsl_vector *yvecmet=gsl_vector_calloc(Rnum+S);

    
  for(l=0;l<Y;l++)						// start of patch solving
  {
    /* Initialisierungen */
    gsl_vector_set_zero(svec);
    

    /* yfdvec enthält die Population */
    gsl_vector_view y_vec=gsl_vector_subvector(yfdvec,(Rnum+S)*l,(Rnum+S));
    gsl_vector *yvecm=&y_vec.vector;
    
    gsl_vector_memcpy(yvecmet,yvecm);

    gsl_vector_memcpy(svec,Mvec);
    //printf("svec vorher: %f\n",gsl_vector_get(svec,3));
    gsl_vector_scale(svec,alpha);				// s(i)=alpha*masse^(-0.25) [svec=Respiration bzw. Mortalitaet]
    //printf("svec nachher: %f\n",gsl_vector_get(svec,3));
    gsl_vector_set(yvecmet,0,0);				// Resource kann nicht aussterben
    gsl_vector_mul(svec,yvecmet);				// s(i) = alpha*masse^(-0.25)*y(i)
   
  
   
    metLoss[l] = gsl_blas_dasum(svec);
    //printf("metloss %f\n",metLoss[0]);
  }
  /* Speicher befreien */
  gsl_vector_free(svec);
  
  return 0;
}


int kastenFunktion(double x, double y, double z)
{
  int flag = 0;
  
  if(x>=y&&x<=z)
  {
    flag = 1;
  }
  return(flag);
}

double* functionalDiversity(struct foodweb nicheweb, const double y[], double center[], double intervall[], double* funDiv)
{

  int S 	 	= nicheweb.S;
  int Y 	     	= nicheweb.Y;
  int Rnum 		= nicheweb.Rnum;
 
  
  
  int i,j,l,max;
  double L=10;
  double N = 1000;
  double h= L/N;
  double s;
  double M_prey;
  double C,I;

  for(l=0; l<Y;l++)
  {
//     int k;
//     k=0;
    for(j=0;j<N;j++)
    {
      s=(j)*h-(2);
      M_prey=s;
      max = 0;
      //printf("M_prey : %f\n", M_prey);
      for(i=(Rnum+S)*l+Rnum; i<(Rnum+S)*l+(Rnum+S);i++)
      {	  

	if(y[i]>1e-5)
	{

	  C = log10(center[i]);
	  I = log10(intervall[i]);


	  if(I>=0)
	  {
// 	    if(M_prey>C-I && M_prey< C+I)
// 	    {
// 	      //printf("k : %i\n",k);
// 	      
// 	      k++;
// 	    }

	    max += kastenFunktion(M_prey,C-I,C+I);
// 	    printf("In Schleife für positive Intervalle\n");
//  	    printf("kastenFunktion : %i\n", kastenFunktion(M_prey,C+I,C-I));
//  	    printf("max : %i\n", max);
	  }
	  else
	  {
// 	    if(M_prey<=C-I && M_prey>= C+I)
// 	    {
// 	      //printf("k : %i\n",k);
// 	      k++;
// 	    }
	    max += kastenFunktion(M_prey,C+I,C-I);
// 	    printf("i : %i\n", i);
// 	    if(kastenFunktion(M_prey,C+I,C-I)>0)
// 	    {
	    
// // 	    printf("M_prey : %f\n", M_prey);
// 	    }
	  }   	    
	 }
	
      }
      //printf("k : %i\n",k);
      //printf("max außen : %i\n", max);
      if(max>0)
      {
	max =1;
      }
      //printf("max Ende: %i\n", max);
      funDiv[l] += max*h;
      //printf("funDiv: %f\n",funDiv[0]);
      }
      
    }
  
  return 0;
}





