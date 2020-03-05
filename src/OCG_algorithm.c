/*************************************************************************/
/*         Program name: Overlapping Class Generator, OCG                */
/*                                                                       */
/*                   Copyright (C) 2011 Alain Guenoche                   */
/*                         guenoche@iml.univ-mrs.fr                      */
/*                                                                       */
/*                        with contributions from:                       */
/*                Benoit Robisson and Charles E. Chapple                 */
/*                                                                       */
/*                                                                       */
/*                                                                       */
/* This program is free software: you can redistribute it and/or modify	 */
/* it under the terms of the GNU General Public License as published by	 */
/* the Free Software Foundation, either version 3 of the License, or	 */
/* (at your option) any later version.					 */
/*   									 */
/* This program is distributed in the hope that it will be useful,	 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of	 */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	 */
/* GNU General Public License for more details.				 */
/*   									 */
/* You should have received a copy of the GNU General Public License	 */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                       */
/*************************************************************************/

/***************************************************************************************************/
/* This program builds an overlapping class system from an unweighted simple graph G=(V,E). 	   */
/* Let |V|=n and |E|=m.										   */
/*  												   */
/* It is essentially a hierarchical ascending algorithm joining two classes at each step.	   */
/* The optimized criterion is the modularity. It can be either the average gain or the total gain. */
/*  												   */
/* The initial overlapping class system can be : 						   */
/* - the set of all maximal cliques (it can take a long time to establish)       		   */
/* - the set of edges (many initial classes (m) implying many steps (O(m))			   */
/* - the set of "centered cliques" (at most n), giving a fast solution for large graphs.	   */
/*  												   */
/* Two class systems can be calculated, the one maximazing the modularity, or the final one.	   */
/* In that case, the expected minimum number of clusters and the maximum caldinality of the 	   */
/* final clusters are required.									   */
/*  												   */
/* Fusion of classes are realized until one of these conditions is fullfiled. 			   */
/* When no more class fusion can be realized the algorithm stops.				   */
/*  												   */
/* Let p be the number of initial classes. The complexity is this algorithm is O(p^3).		   */
/*  												   */
/***************************************************************************************************/
 


/**************************************************** 
 * Code made R-compatible by Alex T. Kalinka, 2012. *
 ****************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <stdint.h>
#include <inttypes.h>


#include <R.h>


/*
  This program computes a hierarchy of overlapping classes 
  The initial classes can be either the maximal cliques of the graph (i), the edges of the graph (ii) or the centred cliques (iii)
  The fusion of two classes optimizing the modularity of teh resulting class system is performed at each step
*/

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/



// Count the number of lines in a file
int CountLines(FILE *fp) {
  int nl = 0;
  char buf[BUFSIZ];
  while ( fgets(buf, sizeof buf, fp) )
    if ( strchr(buf, '\n') ) ++nl;
  if ( !strchr(buf, '\n') ) ++nl;
  return nl;
}



static int compar(const void *e1, const void *e2){
	return strcmp(e1,e2); 
}



int Inclus(short **Cl, short *Q, int N, int *NbCliq)
/*	To search if array Q is included in one of the cliques in Table Cl*/
/*	Sinon on l'ajoute */

{	int k, j, flag;
	
  for (k=0; k < *NbCliq; k++){
	if (Cl[k][0]<0) continue;
      flag=1;
      for (j=0; j < N; j++)
	if (Cl[k][j]<Q[j]) { flag=0; break; }
      if (flag==1) return 1;
    }	// Non inclus
  for (j=0; j < N; j++)
    Cl[*NbCliq][j]=Q[j];
  *NbCliq = *NbCliq + 1;
  return 0;
}



long long Modularity(int **B, int N, char **A)
/*	Computing the integer modularity value */
{	
  int i, j;
  long long Sum=0;

  for (i=1; i<N; i++)
    for (j=0; j<i; j++)
      Sum += A[j][i]*B[i][j];
  return Sum;
}



int compDoubles(double a, double b)
	{
	// To compare doubles without rounding issues.
	double diff;

	diff = a-b;

	if(diff > 0.0000000000000001){
		return 1;
	}else{
		return 0;
		}
	}



void getOCGclusters(char **file, int *ICS, int *FM, int *MCC, int *CCS, int *MC, int *numnodes, int *verb, int *success){


  int   i,j,k, ii, jj = 0, kk = 0, flag, cl1, cl2, k1, k2, ko, mycomp; 
  int 	VarMod, fus1, fus2, card=1;
  long long	BestMod=0;
  double	MMod, var, VarMax;

  int_fast32_t Mod;


  char FichE[120], **A;
  double SumDg2 = 0.0;
  short *Dg, **Cl, *Clas, *Kard, *Q, **BestCl, *BestKard;
  int **B, N=0, Na=0, NbClIni=0, NbClas=0, DgMax=0, BestNbClas=0, CardMax, typ, FuStyl, FCS, ClCh, verbose;
  long **Var;
  long long ModMax=0;

  #define MaxProt 10000 // Maximum number of vertices (Proteins)
  #define SupCar 21

  char Et[*numnodes][SupCar];

  // Set up user-defined variables.
  strcpy(FichE, *file);
  typ = *ICS;
  FuStyl = *FM;
  CardMax = *MCC;
  FCS = *CCS;
  ClCh = *MC;
  verbose = *verb;

   /*************************************************/
   /* Class System should be the Final Class System */
   /* unless centered cliques were chosen as the    */
   /* Initial Class System  			    */
   /*************************************************/
  if(FCS == 1 && typ!=3){
    Rprintf("Argument error!\n\n"); 
    Rprintf("Choice of class system should be 'Final Class System' unless \nCentered Cliques were chosen as the Initial Class System\n\n");
    return;
  }

  /* Setup Errors : Incompatible options selected */

  if (typ!=1 && typ!=2 && typ!=3){
    Rprintf("Argument error!\n\n");
    Rprintf("Type of initial classes: maximal cliques (1) edges (2) centered cliques (3), default = 3\n\n");
    return;
  }
  if (FCS==0){ 
    if(CardMax==0 && ClCh==0){
      Rprintf("Argument error!\n\n"); 
      Rprintf("The Final Class System requires either limiting the maximum class cardinality \n");
      Rprintf("or, setting a minimum number of expected classes\n\n");
      return;
    }
  }
  if (FCS!=0 && FCS!=1){
    Rprintf("Argument error!\n\n"); 
    Rprintf("Final Class System (0) or the one maximizing the modularity (1), default = 0\n\n");
    return;
  }
  
  if (FuStyl!=0 && FuStyl!=1){
    Rprintf("Argument error!\n\n"); 
    Rprintf("Fusion Style must be either 0 (according to the average) or 1 (total gain), default = 0\n\n");
    return;
  }

  if(typ==1 || typ==2){
   if(FCS==1){
    Rprintf("Argument error!\n\n"); 
    Rprintf("The final Class System maximizing the modularity can only be used with the centered cliques\n\n");
    return;
    }
  }

  if(FCS==1 && CardMax!=0){
    Rprintf("Argument error!\n\n"); 
    Rprintf("The final Class System maximizing the modularity can only be used without setting the maximum cardinality for each class\n\n" );
    return;
  }
    
  if(FCS==1 && ClCh!=2){
    Rprintf("Argument error!\n\n"); 
    Rprintf("The final Class System maximizing the modularity can only be used with a minimum number of expected classes of 2\n\n");
    return;
  }


  /*  Establish a hierarchy of overlapping classes optimizing the modularity criterion  */


  int	NR;
  char	Ch1[2*SupCar],Ch2[2*SupCar],OldCh[SupCar]="";

  char buf[BUFSIZ];
  int NbCar;
  char MaxCar=0;

  FILE *FichCar;

  /********************************************************************/
  /* Read an unweighted graph file, defined by 			      */
  /* the number of records (NR) and NR records containing	      */
  /* Each edge is defined by two labels (corresponding to 	      */
  /* vertices) ; the same edge can be duplicated with no effect	      */
  /* Labels are limited to 20 characters			      */
  /********************************************************************/


  FichCar = fopen(FichE,"r");
  if ( FichCar == 0){
    Rprintf( "Could not open file: %s\n", FichE);
    return;
  }
 

  // Read in file.

  NR=CountLines(FichCar); 

  rewind(FichCar); 
  fclose(FichCar);

  FichCar = fopen(FichE,"r");

  strcpy(OldCh,""); NbCar=0;
  for (i=0; i<NR; i++){

	if( fgets(buf, sizeof buf, FichCar) != NULL){
		sscanf(buf, "%s %s", Ch1, Ch2);
	}else{
		return;
		}

      if (strlen(Ch1)>NbCar) NbCar=strlen(Ch1); 
      if (strlen(Ch2)>NbCar) NbCar=strlen(Ch2); 
      if (NbCar>MaxCar) MaxCar=NbCar;
      if (MaxCar>=SupCar){
	Rprintf ("Labels are limited to %d characters\n",(SupCar-1));
	  Rprintf("Incorrect edge: %5d   %s, %s\n",i,Ch1,Ch2);  
	  return;
        }
      if (strcmp(Ch1,OldCh) != 0) // if it is not the same as before.
    	{	flag=1;
	  for (j=0; j < N; j++)
	    if (strcmp(Ch1,Et[j]) == 0) { flag=0; strcpy(OldCh,Ch1); break; }
	  if (flag){
	strcpy(Et[N],Ch1); N++; 
	      if (N==MaxProt) { Rprintf("warning: more than %d vertices; algorithm will take > 1 hour to run",MaxProt);}
	    }	}
      flag=1;
      for (j=0; j < N; j++)
	if (strcmp(Ch2,Et[j]) == 0) { flag=0; break; }
      if (flag){
	strcpy(Et[N],Ch2); N++; 
	  if (N==MaxProt) {  Rprintf("warning: more than %d vertices; algorithm will take > 1 hour to run ",MaxProt);}
    	}

    }	

  rewind(FichCar);
  fclose(FichCar);

  qsort(Et, N, SupCar, compar);

  /*	Second lecture to store the graph  */

  A = malloc(N * sizeof(char *));

  assert (A != NULL);
  for (i=0; i < N; i++){
      A[i] = malloc ( N * sizeof(char) );	/* A[i][j]=1 if i adjacent a j, A[i][j]=0 otherwise */
      assert ( A[i] != NULL);
      for (j=0; j < N; j++) A[i][j]=0;
    }
  /*	A is the two dimensional adjacency table : A[i][j]=A[j][i] iff vertices i and j are linked 
	B is the two dimensional table containing the (integer) modularity value of any pair (i,j)
	Clas is a one dimensional array : Clas[i] is the class number of vertex i
	Dg is a one dimensional array : Dg[i] is the degree (nb. of adjacent vertices) of vertex i
  */	

  FichCar = fopen(FichE,"r");  

 if ( FichCar == 0){
   Rprintf( "Could not open file: %s\n", FichE);
    return;
  }
 NR=CountLines(FichCar);

  rewind(FichCar);
  fclose(FichCar);
  
  FichCar = fopen(FichE,"r");  

 if ( FichCar == 0){
    Rprintf( "Could not open file: %s\n", FichE );
    return;
  }
  //  fscanf(FichCar,"%d",&NR);
  strcpy(OldCh,""); ii=-1;
  for (i=0; i<NR; i++){ 

	if( fgets(buf, sizeof buf, FichCar) != NULL){
		sscanf(buf, "%s %s", Ch1, Ch2);
	}else{
		return;
		}

      if (strcmp(Ch1,Ch2) == 0) continue;
      if (strcmp(Ch1,OldCh) != 0) // if it is not the same as before
    	{	for (j=0; j < N; j++)
	    if (strcmp(Ch1,Et[j]) == 0) { ii=j; strcpy(OldCh,Ch1); break; }
    	}
      for (j=0; j < N; j++)
	if (strcmp(Ch2,Et[j]) == 0) { jj=j; break; }
      A[ii][jj]=1; A[jj][ii]=1;

    }
  rewind(FichCar);
  fclose(FichCar); 


  B = malloc( N * sizeof(int *));
  Dg = malloc( N * sizeof(short));

  for(i=0; i < N; i++){
	B[i] = malloc ( N * sizeof(int) );
	assert ( B[i] != NULL);
	}
//Initialise.
for(i=0;i<N;i++){
 Dg[i]=0;
 for(j=0;j<N;j++){
  B[i][j]=0;
  }
 }


  /* Evaluating the number of edges (Na) the vertex degrees, and DgMax its maximum value*/
  DgMax=0; Na=0; 
  for (i=0; i < N; i++){
   Dg[i]=0;
      for (j=0; j < N; j++)
	if (A[i][j] > 0) Dg[i] = Dg[i] + 1;
      if (Dg[i] > DgMax) DgMax = Dg[i];
      Na += Dg[i]; // it is the double
      SumDg2 = SumDg2 + (double) Dg[i]* (double) Dg[i];
    }
//Rprintf("SumDg2 = %5f\n",SumDg2);


  /* Evaluating the Modularity values and its maximum */
  ModMax=0;
  for (i=1; i < N; i++)
    for (j=0; j<i; j++){
	B[i][j] = Na*A[i][j]-Dg[i]*Dg[j];
	B[j][i] = B[i][j];
	ModMax += A[i][j]*B[i][j];
      }


  Na = Na/2;

  Clas = malloc(N * sizeof(short));

//Initialise.
for(i=0;i<N;i++)
 Clas[i]=0;


    if (FCS==1){
      CardMax=N;
      ClCh=2;
    }
    else 
      {
	if (CardMax==0) CardMax=N;
	if (ClCh==0) ClCh=2;
      }

  int Adj, NbAdj, mis=0, Somcard=0;
  short *D;

  switch(typ){

    case 1:
	Cl = malloc((2*N) * sizeof(short *));
	Kard = malloc(N * sizeof(short));
	Q = malloc(N * sizeof(short));
      //NbClIni = Clique(Cl,Kard,Q,N,A);
  int   NbY, Nl, Long, NewLong, cli, flag1, flag2, flag3;
  short	*X, *Y;
  int NbCliq = 1;

//  Q = malloc((N) * sizeof(short));
  X = malloc(N * sizeof(short));
  Y = malloc(N * sizeof(short));
  assert(Q != NULL && X != NULL && Y != NULL);


  /*  Cl is a two dimensional table that will contain all the maximal cliques as initial classes
      Each row corresponds to a clique
  */

  Nl=2*N; Long=Nl; // Allocating Nl rows
//  Cl = malloc((Nl) * sizeof(short *));
  assert(Cl != NULL);
  for (i=0; i<Nl; i++){
	Cl[i] = malloc (N * sizeof(short));
        assert (Cl[i] != NULL);
    }
// Initialise.
for(i=0;i<N;i++){
	Q[i]=0;X[i]=0;Y[i]=0;Kard[i]=0;
}
for(i=0;i<Nl;i++)
 for(j=0;j<N;j++)
	Cl[i][j]=0;


  Cl[0][0]=1;
  for (j=1; j<N; j++)
    if (A[j][0]) Cl[0][j]=1; else Cl[0][j]=0;
  NbCliq=1; 
	
  //	Sequential algorithm described in the book
  for (i=1; i<N-1; i++){
	for (j=0; j<i; j++) X[j]=0;
      X[i]=1;
      for (j=i+1; j<N; j++)
	if (A[j][i]) X[j]=1; else X[j]=0;
		
      // Allocating again
      cli=0;
      for (k=0; k < NbCliq; k++)
	if (Cl[k][i]>0) cli++;
      if (NbCliq + 2*cli+1 > Long){
	  NewLong=Long+Nl;
	  Cl = realloc(Cl,NewLong*sizeof(short *));
	  assert (Cl != NULL);
	  for (k=Long; k<NewLong; k++){
		Cl[k] = malloc(N * sizeof(short));
		assert (Cl[k] != NULL);
	    }
	  //Initialise.
	  for(j=Long;j<NewLong;j++)
	   for(k=0;k<N;k++)
	    Cl[j][k]=0;
	Long = NewLong;
	}
      for (k=0; k < NbCliq; k++){
	if (Cl[k][i]==0) continue;
	  NbY=0;
	  for (j=i+1; j<N; j++)
	    if (Cl[k][j]==1 && X[j]==0) { Y[NbY]=j; NbY++; }
	  if (NbY>0){
	for (j=0; j<N; j++)
		Q[j]=Cl[k][j];
	      Cl[k][0]=-1;
	      Q[i]=0; flag1 = Inclus(Cl,Q,N,&NbCliq);
	      Q[i]=1;
	      for (j=0; j<NbY; j++)
		Q[Y[j]]=0;
	      flag2 = Inclus(Cl,Q,N,&NbCliq);
	      if (flag1+flag2 == 2) continue;
	      for (j=0; j<N; j++) 
		Cl[k][j]=Cl[NbCliq-1][j];
	      NbCliq--;
	    }	
	}
      for (j=0; j<N; j++) Q[j] = X[j];
      flag3 = Inclus(Cl,Q,N,&NbCliq);
    }
  free(Kard);
  Kard = malloc(NbCliq * sizeof(short));
  //Initialise.
  for(i=0;i<NbCliq;i++)
	Kard[i]=0;
//  assert(Kard != NULL);
  //	Eliminating the unused rows and coding cliques as item listes
  k=0;
  for (i=0; i < NbCliq; i++){
   if (Cl[i][0]<0) continue; 
      kk=0;
      for (j=0; j<N; j++)
	if (Cl[i][j]) { Cl[k][kk]=j; kk++; }
      Kard[k]=kk;
      k++;
    }  NbCliq=k;
  free(Q); free(X); free(Y);
  NbClIni = NbCliq;

      break;

    case 2:
	Cl = malloc(Na * sizeof(short *));
	Kard = malloc(Na * sizeof(short));
      //NbClIni = ClasArete(CardMax,Cl,Kard,N,Na,A);

  /*  The initial class system is the set of edges
      Cl is a two dimensional table that will contain all the initial classes
      Kard is an array containing the number of elements in each class
  */
//  Cl = malloc((Na) * sizeof(short *));
//  Kard = malloc((Na) * sizeof(short));
  assert(Cl != NULL && Kard != NULL);
  for (i=0; i<Na; i++){
	Cl[i] = malloc (CardMax * sizeof(short));
	assert (Cl[i] != NULL);
	Kard[i]=2;
    }
// Initialise.
for(i=0;i<Na;i++)
 for(j=0;j<CardMax;j++)
  Cl[i][j] = 0;

  k=0;
  for (i=1; i < N; i++)
    for (j=0; j<i; j++)
      if (A[i][j]){
   Cl[k][0]=j; Cl[k][1]=i; k++; }
   NbClIni = k;

      break;

    case 3:

  /*  Compute the centered cliques as initial class system */
  /*  Cl is a two dimensional table that will contain all the centered cliques as initial classes
      Kard is an array containing the number of elements in each clique
      D is a one dimentional array indicating the internal degre of each vertex (nb. of adjacent vertices in its class)
      X is a working vector
  */
  Cl = malloc(N * sizeof(short *));
  Kard = malloc(N * sizeof(short));
  X = malloc(N * sizeof(short));
  D = malloc(N * sizeof(short)); 
  assert(Cl != NULL && Kard != NULL && X != NULL && D != NULL);
  for (i=0; i < N; i++){
	Cl[i] = malloc (N * sizeof(short));
	assert (Cl[i] != NULL);
	X[i]=0;
    }
//Initialise.
for(i=0;i<N;i++){
 Kard[i] = 0;
 for(j=0;j<N;j++){
  Cl[i][j]=0;
 }
}


  if(verbose==1){Rprintf("Calculating Initial class System...");}
  for (i=0; i < N; i++){

      if(verbose==1 && i%200==0){Rprintf(".");}
      
      Cl[NbClas][0]=i; 
      card=1; NbAdj=0;
      for (j=0; j < N; j++)
	if (A[i][j]==1) { X[NbAdj]=j; D[NbAdj]=1; NbAdj++; }
      // Calcul des degres relatifs a X 
      for (k=1; k<NbAdj; k++)
	for (j=0; j<k; j++)
	  if (A[X[k]][X[j]]==1) { D[k]++; D[j]++; }

      for (Adj=0; Adj<NbAdj; Adj++){
	DgMax=0;
	  for (k=0; k<NbAdj; k++)
	    if (D[k] >= DgMax) { DgMax=D[k]; kk=k; }
	  if (DgMax < card) break; // Not enough adjacent vertices
	  jj=X[kk]; D[kk]=0;
	  kk=0;  
	  for (k=0; k<card; k++){
		j=Cl[NbClas][k];
	      if (A[jj][j]==1) kk++;
	    }
	  if (kk-card >= 0) { Cl[NbClas][card]=jj; card++; } 
	}
      Kard[NbClas]=card;
      // Is it a new class ?
      flag=0;
      for (k=0; k<NbClas; k++){
	if (Kard[k]<card) continue;
	  for (j=0; j < N; j++) X[j]=0;
	  for (j=0; j<Kard[k]; j++) X[Cl[k][j]]=1;
	  flag=1;
	  for (j=0; j<Kard[NbClas]; j++) if (X[Cl[NbClas][j]]==0) flag=0;
	  if (flag==1) break;  // The new class exists
	}
      if (flag==0){
	Somcard += card;
	NbClas++;
	}
    }


  if(verbose==1){Rprintf("Done\n");}

  free(X); free(D);

  if(verbose==1){Rprintf("Nb. of classes %d\n",NbClas);}

  //mis = Fermeture(NbClas,Cl,Kard,N,A);

  /* This procedure updates the upper right part of Table A
     When i > j, A[i][j]=1 iff vertices i and j are joined in at least one class
     When i < j, A[i][j]=1 iff there is an edge between i and j
  */
  for (i=1; i < N; i++)
    for (j=0; j<i; j++)
      A[j][i]=0;	
  for (i=0; i<NbClas; i++)
    for (j=1; j<Kard[i]; j++){
   jj=Cl[i][j];
	for (k=0; k<j; k++){
   kk=Cl[i][k];
	    if (kk>jj) A[jj][kk]=1; else if (kk<jj) A[kk][jj]=1;
	  }
      }


  // Count the edges such the two ends are not joined in any class
  for (i=1; i < N; i++)
    for (j=0; j<i; j++)
      if (A[i][j]==1 && A[j][i]==0) mis++;

  if(verbose==1){Rprintf("Nb. of edges not within the classes %d\n",mis);}

  NbClIni = NbClas;

    }

  if(verbose==1){Rprintf("Number of initial classes %d\n",NbClIni);}



  //Mod = Modularity(B,N,A); 

  Mod=0;

  for (i=1; i<N; i++)
    for (j=0; j<i; j++)
      Mod += A[j][i]*B[i][j];


  MMod = 1.*Mod/Na/Na/2-SumDg2/Na/Na/4;

//Rprintf("Mod = %d, MMod = %5f\n",Mod,MMod);

//  if (N<20) ClasOut(NbClIni,Mod,MMod,Cl,Clas,Kard,N,Na,FuStyl,typ,CardMax,ClCh,FCS,FichE,A);


  /**********************************************************************************************/
  /* Matrix Var contains two types of information, one in the upper right and the other	        */
  /* in the bottom left. 								        */
  /* 											        */
  /* When i > j, Var[i][j] is the modularity variation given by the fusion of classes  i and j. */
  /* When i < j, Var[i][j] = 1 if they can be merged and 0 if not.			        */
  /* Two classes can be merged if they are connected and if their fusion does not pass CardMax. */
  /**********************************************************************************************/

  Var = malloc((NbClIni) * sizeof(long *));
  BestCl = malloc((NbClIni) * sizeof(short *));
  BestKard = malloc((NbClIni+1) * sizeof(short));
  assert (Var != NULL && BestCl != NULL && BestKard != NULL);
  if(verbose==1){Rprintf("Running...");}
  for (i=0; i<NbClIni; i++){	
	Var[i] = malloc ( NbClIni * sizeof(long) );	
	BestCl[i] = malloc ( CardMax * sizeof(short) );
	assert ( Var[i] != NULL && BestCl[i] != NULL);
	Var[i][i]=-1;
      for (j=0; j<i; j++) { Var[i][j]=-Mod; Var[j][i]=0; }
    }

  //	Initialization of the Var table comparing any two initial classes
  for (cl1=1; cl1<NbClIni; cl1++){

      if(verbose==1 && cl1%100==0){Rprintf(".");}
      for (cl2=0; cl2<cl1; cl2++){
   	VarMod=0; fus1=0;
	  for (i=0; i<N; i++) Clas[i]=0;
	  for (k1=0; k1<Kard[cl1]; k1++){
   		i=Cl[cl1][k1]; Clas[i]=1;
	      for (k2=0; k2<Kard[cl2]; k2++){
  		 j=Cl[cl2][k2]; Clas[j]=1;
		  if (i==j) fus1=1;  // Classes connexes
		  else if (i<j && A[j][i]==1) fus1=1;  // Connected classes 
		  else if (i>j && A[i][j]==1) fus1=1;  // Connected classes
		  if (i<j) VarMod += (1-A[i][j])*B[i][j]; else if (i>j) VarMod += (1-A[j][i])*B[i][j]; 
		}
	      kk=0;
	      for (k=0; k<N; k++) if (Clas[k]) kk++;
	      if (kk > CardMax) fus1=0;
	    } 
	  Var[cl1][cl2] = VarMod; Var[cl2][cl1] = fus1;
	}
    }
  if(verbose==1){Rprintf(".\n");}
  NbClas = NbClIni; card=1;




  while (NbClas > ClCh){

      VarMax = -1*(double)ModMax; cl1=-1;
      for (ii=1; ii<NbClIni; ii++){
   	if (Kard[ii]==0) continue;		// old classes are indicated by setting their cardinality to 0
	  for (jj=0; jj<ii; jj++){
   		if (Kard[jj]==0) continue;
	      // Can we fuse them?
	      if (Var[jj][ii]==0) continue;		// no 
	      if (FuStyl==0) card = Kard[ii]*Kard[jj]; 	// average gain
	      var = (double)Var[ii][jj]/(double)card;
	      mycomp = compDoubles(var,VarMax);
	      if (mycomp == 1) { VarMax = var; cl1=jj; cl2=ii; }
	    }

	}

      if (cl1<0) {  break; }
      for (k1=0; k1 < Kard[cl1]; k1++){
	i=Cl[cl1][k1]; 
	  for (k2=0; k2 < Kard[cl2]; k2++){
		j = Cl[cl2][k2];
	      // The joined pairs are indicated in the lower left part of A
	      if (i<j) A[i][j]=1; else if (i>j) A[j][i]=1;
	    }
	}
      //Mod = Modularity(B,N,A);

  Mod=0;

  for (i=1; i<N; i++)
    for (j=0; j<i; j++)
      Mod += A[j][i]*B[i][j];


      for (k=0; k<N; k++) Clas[k]=0;
      for (k1=0; k1<Kard[cl1]; k1++) Clas[Cl[cl1][k1]]=1;
      for (k2=0; k2<Kard[cl2]; k2++) Clas[Cl[cl2][k2]]=1;
      kk=0;
      for (k=0; k<N; k++) if (Clas[k]) { Cl[cl1][kk]=k; kk++; }
      Kard[cl1]=kk; Kard[cl2]=0;  // Cl[cl1][] will contain the new class
  
      if (Mod >= BestMod){
   		kk=0;
	  for (k=0; k<NbClIni; k++){
   		if (Kard[k]==0) continue;
	      for (i=0; i<Kard[k]; i++)
		BestCl[kk][i]=Cl[k][i];
	      BestKard[kk]=Kard[k];
	      kk++;
	    }
	  BestNbClas = NbClas-1;
	  BestMod = Mod;
	}
      // Updating table Var
      for (ii=0; ii<NbClIni; ii++){
	if (Kard[ii]==0 || cl1==ii || cl2==ii) continue;
	  fus1=0; fus2=0;
	  if (ii < cl1) fus1=Var[ii][cl1]; else fus1=Var[cl1][ii];
	  if (ii < cl2) fus2=Var[ii][cl2]; else fus2=Var[cl2][ii];
	  if (fus1==0 && fus2==0) // Can be fused with neither one nor the other
	    {   if (ii < cl1) Var[ii][cl1]=0; else Var[cl1][ii]=0;
	      continue;
	    }
	  // verifying the cardinality condition and computing Var[ii][cl1]
	  VarMod=0; 
	  for (k=0; k<N; k++) Clas[k]=0;
	  for (k=0; k<Kard[ii]; k++){
	    i=Cl[ii][k]; Clas[i]=1;
	      for (k1=0; k1<Kard[cl1]; k1++){
		j=Cl[cl1][k1]; Clas[j]=1;
		if (i<j){ VarMod += (1-A[i][j])*B[i][j];}
		else if (i>j){ VarMod += (1-A[j][i])*B[i][j];}
	      }
	  }   
	  kk=0;	// cardinality of the eventual fusion
	  for (k=0; k<N; k++) if (Clas[k]) kk++;
	  if (kk > CardMax) fus1=0; else fus1=1;
			
	  if (ii > cl1) { Var[ii][cl1]=VarMod; Var[cl1][ii]=fus1; }
	  else		{ Var[cl1][ii]=VarMod; Var[ii][cl1]=fus1; }
	}
      NbClas--; 
      if (NbClas%10 == 0 && verbose==1) {Rprintf("Remaining classes: %d of %d     \r",NbClas,NbClIni);
					 R_FlushConsole();
					 R_ProcessEvents();}
    }


  
  Rprintf("Remaining classes: None                 \n");

  if (FCS==1)
    //	Going back to classes making the largest modularity
    {
      for (k=0; k<BestNbClas; k++){
	Kard[k] = BestKard[k];
	  for (i=0; i<Kard[k]; i++)
	    Cl[k][i] = BestCl[k][i]; // printf("%s ",Et[Cl[k][i]]);
	}   NbClas = BestNbClas;
      for (k=NbClas+1; k<NbClIni; k++) Kard[k]=0;
    }
  else
    {   kk=0;
      for (k=0; k<NbClIni; k++){
   if (Kard[k]==0) continue;
	  for (i=0; i<Kard[k]; i++){
	    Cl[kk][i] = Cl[k][i];}
	  Kard[kk] = Kard[k];
	  kk++;
	}
      NbClas=kk;
    }



  //  Making singletons with vertices having a negative contribution to the mudularity of their class 

  //Effacer(Cl,Clas,Kard,B,N,NbClas);

  int		ncard, NbElEl=0;
  long long	Cont;


  
  // Calculate the contribution of each element within its class
  for (k=0; k<NbClas; k++){
   
      for (i=0; i<N; i++) Clas[i]=0;
      for (i=0; i<Kard[k]; i++){
	  ii = Cl[k][i]; //printf("%d ",ii+1);
	  Cont=0;
	  for (j=0; j<Kard[k]; j++){
		jj = Cl[k][j];
		Cont += B[ii][jj];
	    }
	  if (Cont >= 0) Clas[ii]=1; 
	  else { 
	    NbElEl++; }
	}
      ncard=0;
      for (i=0; i<N; i++)
	if (Clas[i] == 1) { Cl[k][ncard] = i; ncard++; }
      Kard[k] = ncard; 
		
    }

//   int ko=0;
//   for (i=1; i<N; i++){
//     for (j=0; j<i; j++){
//       if (A[i][j] && A[j][i]){ ko++;}
//     }
//   }
//Rprintf("NbClas = %d\nko = %d\n",NbClas,ko);



  for (i=0; i<N-1; i++)
    for (j=i+1; j<N; j++)
      A[i][j]=0;
  for (k=0; k < NbClas; k++)
    for (i=0; i < Kard[k]-1; i++){
	ii = Cl[k][i];
	for (j=i+1; j < Kard[k]; j++){
	jj = Cl[k][j];
	    if (ii < jj){ A[ii][jj]=1; }else{ A[jj][ii]=1;}
	  }
      }

//   ko=0;
//   for (i=1; i<N; i++){
//     for (j=0; j<i; j++){
//       if (A[i][j] && A[j][i]){ ko++;}
//     }
//   }
//Rprintf("NbClas = %d\nko = %d\n",NbClas,ko);
   //Mod = Modularity(B,N,A); 

  Mod=0;

  for (i=1; i<N; i++)
    for (j=0; j<i; j++)
      Mod += A[j][i]*B[i][j];


  MMod=1.*Mod/Na/Na/2 - 1.*SumDg2/Na/Na/4;


   FILE *OUT;
   OUT = fopen("OCG_temp.txt","w");


   fprintf(OUT, "#################################################\n");
   //   printf("# %s results on graph file %s\n", argv[0],FichE);
   fprintf(OUT, "# Graph file %s has %d Vertices and %d edges\n",FichE,N,Na);
   //fprintf(OUT, "# Degree maximum: %d\n",DgMax);
   fprintf(OUT, "# Rate of edges: %.4f\n",1.*Na/N/(N-1));
   
   ko=0;
   for (i=1; i < N; i++){
     for (j=0; j<i; j++){
       if (A[i][j] && A[j][i]){ ko++;}
     }
   }
   fprintf(OUT, "# Rate of intraclass edges: %.3f\n#\n",1.*ko/Na);

   
 switch(typ){
    case 1:
      fprintf(OUT, "# Initial classes are: maximal cliques\n");
      break;
    case 2:
      fprintf(OUT, "# Initial classes are: edges\n");
      break;
    case 3:
      fprintf(OUT, "# Initial classes are: centered cliques\n");
    }	


   fprintf(OUT, "# Fusion is according to: ");
   if (FuStyl==0) { fprintf(OUT, "the average gain\n");}
   else{ fprintf(OUT, "the total gain\n");}
   if (FCS==1){
     fprintf(OUT, "# Class System is: maximizing the modularity\n");
   }
   else{
     fprintf(OUT, "# Class System is: Final Class System\n\n");
     fprintf(OUT, "# Maximum Class Cardinality: %d\n",CardMax);
     fprintf(OUT, "# Minimum Number of Classes %d\n",ClCh);
   }
   fprintf(OUT, "#################################################\n");

  for (i=0; i < N; i++) Clas[i]=0;
  for (k=0; k<NbClas; k++){
      for (i=0; i<Kard[k]; i++){
	  Clas[Cl[k][i]]++;
	}
    }
  k=0; kk=0;
  for (i=0; i < N; i++) if (Clas[i]==0) kk++; else if (Clas[i]>1) k++; 
  if (kk){
    fprintf(OUT, "\nUnclustered nodes (%d): \n",kk);
    for (i=0; i < N; i++) 
	if (Clas[i]==0) fprintf(OUT, "%s  ",Et[i]);
      fprintf(OUT, "\n");
  }
  else{
    fprintf(OUT, "\nUnclustered nodes (%d): None\n",kk);
  }
  fprintf(OUT, "Multiclustered nodes (%d):\n",k);
  for (i=0; i < N; i++){ 
    if (Clas[i]>1) { 
      fprintf(OUT, "%s (%d), ",Et[i],Clas[i]); 
    }
  }
  fprintf(OUT, "\n");

  /* Now print the classes */
  fprintf(OUT, "\n\nFinal classes (%d), Modularity = %"PRIdFAST32" (%.4f):\n", NbClas, Mod, MMod );
  for (k=0; k<NbClas; k++){

      fprintf(OUT, ">Class %3d  (%d nodes):\n",k+1,Kard[k]);

      for (i=0; i<Kard[k]; i++){
  
	  fprintf(OUT, "%s ",Et[Cl[k][i]]);
	} fprintf(OUT, "\n");
    } fprintf(OUT, "\n");
  fclose(OUT);



  // Clean up memory.
  free(Clas); free(Kard);
  free(BestKard); free(Dg);

  for (i = 0; i < N; i++) { 
  free(A[i]);
  free(B[i]);
  free(Cl[i]);
  }
  free(A);
  free(B);
  free(Cl);

  for (i = 0; i < NbClIni; i++) { 
  free(Var[i]);
  free(BestCl[i]);
  }
  free(Var);
  free(BestCl);


  *success = 1;

}




