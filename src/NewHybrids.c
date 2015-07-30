#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "ECA_ReadData.h"
#include "ECA_MemAlloc.h"
#include "MCTypesEtc.h"
#include "MathStatRand.h"
#include "ranlib.h"


#include "NewHybrids.h"

/* 
	Given total counts N[], and number of "1"'s y[], at distance x[] for each of 
	location which HasLocale and for which the largest index is MaxIndex; and for
	logistic parameters a and b: Prob(1|x) = exp(b(x-a)) / (1+exp(b(x-a))), 
	return the log likelihood function.
*/
double LogistLogL(double a, double b, ival **N, ival **y, double *x, int *HasLocale, int MaxIndex)
{
	int i;
	double p,q, expo,ret=0.0;
	
	for(i=0;i<=MaxIndex;i++) {
		if(HasLocale[i]) {
			/* compute p */
			expo = exp(b*(x[i]-a));
			p = expo/(1.0 + expo);
			q = 1.0 - p;
			
			/* put on the terms in the LogL */
			ret += y[i]->v*log(p) + (N[i]->v - y[i]->v)*log(1-p);
			
			/* printf("X=%f   p=%f   q=%f  \n",x[i],p,q); */
		}
	}
	return(ret);
}

/* initialize the cline variables */
/*This basically means choosing alpha and beta from random values */
void RandoStartClineVars(hyb_chain *C, cline_struct *L)
{
	int l;
	locline_struct *Loc;
	int N,lo,hi;
	
	if(L->NumLocales==0) {
		return;
	}
	
	N = L->NumClineXs;
	lo = (int)(.15 * N);
	hi = (int)(.85 * N);
	
	
	CYCLE_l(C->Dat)   /* cycle over loci */
		/* define this pointer to the current locline struct */
		Loc = L->Locs[l];
		

		Loc->alpha->v = L->ClineX[UniformRV(lo,hi)];
		Loc->beta->v = ranf() * .1;
		
		Loc->LogLike->v = LogistLogL(Loc->alpha->v,Loc->beta->v,Loc->N,Loc->Y,L->x,L->HasLocale,L->MaxLocaleIndex);

			
	END1CYCLE
}	

/* here we accumulate the averages of these cline vars */
void IncrementClineVars(hyb_chain *C, cline_struct *L)
{
	int p,l,x;
	
	if(L->NumLocales==0) {
		return;
	}
	
	/* get the averages for Y and N and alpha and beta */
	for(p=0;p<=L->MaxLocaleIndex;p++) {
		if(L->HasLocale[p]) {
			CYCLE_l(C->Dat)
				IncrementIval(L->Locs[l]->Y[p]);
				IncrementIval(L->Locs[l]->N[p]);
				
			END1CYCLE
		}
	}
	
	CYCLE_l(C->Dat)
		IncrementDval(L->Locs[l]->alpha);
		IncrementDval(L->Locs[l]->beta);
		
		for(x=0;x<L->NumClineXs;x++)  {
			IncrementDval(L->Locs[l]->ClineY[x]);
		}

		
	END1CYCLE
}

/* here we reset cline vars averages to zero */
void ResetClineVarAverages(hyb_chain *C, cline_struct *L)
{
	int p,l,x;
	
	if(L->NumLocales==0) {
		return;
	}
	
	/* get the averages for Y and N and alpha and beta */
	for(p=0;p<=L->MaxLocaleIndex;p++) {
		if(L->HasLocale[p]) {
			CYCLE_l(C->Dat)
				InitIvalSummaryToZero(L->Locs[l]->Y[p]);
				InitIvalSummaryToZero(L->Locs[l]->N[p]);
				
			END1CYCLE
		}
	}
	
	CYCLE_l(C->Dat)
		InitDvalSummaryToZero(L->Locs[l]->alpha);
		InitDvalSummaryToZero(L->Locs[l]->beta);
		
		for(x=0;x<L->NumClineXs;x++)  {
			InitDvalSummaryToZero(L->Locs[l]->ClineY[x]);
		}

		
	END1CYCLE
}

/* Cline Variables Update: just a big function to set all the cline parameters. */
void UpdateClineVariables(hyb_chain *C, cline_struct *L)
{
	int i,p,l,a,d;
	int locale,WasUpdated;
	locline_struct *Loc;
	double aprime, bprime, proppi;
	
	/* first we need to count the gene copies in each locale assigned to the "1" species */
	
	/* initialize everything to zero */
	for(p=0;p<=L->MaxLocaleIndex;p++) {
		if(L->HasLocale[p]) {
			CYCLE_l(C->Dat)
				L->Locs[l]->Y[p]->v = 0;
				L->Locs[l]->N[p]->v = 0;
			END1CYCLE
		}
	}
	
	
	/* then cycle over individuals and loci and count up Y and N */
 	CYCLE_i(C->Dat) 
		if(C->Lat->Ind[i]->Flag[HAS_LOCALE]) {
			locale = C->Lat->Ind[i]->Flag[LOCALE];
			
			CYCLE_l(C->Dat)
				if(YnotMissing(C,i,l)) {  /* only to this if Y is not missing */
					if(C->Lat->Ind[i]->Flag[Z_FIXED]) {
						if(C->Lat->Ind[i]->Flag[FIXED_Z]==1)  {  /* note z=1 must be species one for this to work! */
							L->Locs[l]->Y[locale]->v += 2;
						}
						L->Locs[l]->N[locale]->v += 2;
					}
					else {
						CYCLE_a
							if(C->Lat->Ind[i]->W[l][a]->v==1) {
								L->Locs[l]->Y[locale]->v++;
							}
							L->Locs[l]->N[locale]->v++;
						END1CYCLE
					}
				}
			END1CYCLE
		}
	END1CYCLE
	

/* this is for checking/debugging	
	for(p=0;p<=L->MaxLocaleIndex;p++) {
		if(L->HasLocale[p]) {
			CYCLE_l(C->Dat)
			printf("Locale %d  has Y= %d and N=%d at locus %d\n", p,L->Locs[l]->Y[p]->v,L->Locs[l]->N[p]->v,l);
			END1CYCLE
		}
	}
/* */


	/* cool, now we are ready to do the updates for alpha and beta at each locus, so we cycle over loci and do it */
	CYCLE_l(C->Dat)   /* cycle over loci */
		/* define this pointer to the current locline struct */
		Loc = L->Locs[l];
		WasUpdated = 0;
		
		
		/* compute log likelihood for the current set of Y, N, alpha, and beta, and store it */
		Loc->LogLike->v = LogistLogL(Loc->alpha->v,Loc->beta->v,Loc->N,Loc->Y,L->x,L->HasLocale,L->MaxLocaleIndex);

		
		/* propose new a-value */
		aprime = gennor((float)Loc->alpha->v,(float)Loc->a_sd->v);
		if(aprime>0.0) {  /* reject out of hand if less than zero */
			proppi = LogistLogL(aprime,Loc->beta->v,Loc->N,Loc->Y,L->x,L->HasLocale,L->MaxLocaleIndex);
			if(ranf() < exp(proppi-Loc->LogLike->v)) {
				Loc->alpha->v = aprime;
				Loc->LogLike->v = proppi;
				WasUpdated  = 1;
			}
		}
		bprime = gennor((float)Loc->beta->v,(float)Loc->b_sd->v);
		if(bprime>0.0) {  /* reject out of hand if less than zero */
			proppi = LogistLogL(Loc->alpha->v,bprime,Loc->N,Loc->Y,L->x,L->HasLocale,L->MaxLocaleIndex);
			if(ranf() < exp(proppi - Loc->LogLike->v)) {
				Loc->beta->v = bprime;
				Loc->LogLike->v = proppi;
				WasUpdated = 1;
			}
		}
		if(WasUpdated   || 1) {int x;  /* if we updated these things, then we need to recalculate the cline curve itself
											NOTE: I'll just have it update these every time, even if not necessary.  */
			for(x=0;x<L->NumClineXs;x++)  { double xx; double aa; double bb; double TheY;
				xx = L->ClineX[x];
				aa = Loc->alpha->v;
				bb = Loc->beta->v;
				TheY = exp( bb * (xx - aa) ) / ( 1 + exp( bb * (xx - aa) ) );
/*				printf("x = %d  TheX = %f  TheY = %f\n",x,xx, TheY); */
				Loc->ClineY[x]->v = TheY;
			}
		}
	END1CYCLE
	
/* This is just for debugging stuff 
l=1;
printf("alpha = %f   beta = %f\n",L->Locs[l]->alpha->v,L->Locs[l]->beta->v);
for(d=0;d<=L->MaxLocaleIndex;d++) {
	if(L->HasLocale[d]) {
		printf("%d   %d   %.2f    %.1f\n",L->Locs[l]->Y[d]->v,
					L->Locs[l]->N[d]->v,(double)L->Locs[l]->Y[d]->v/(double)L->Locs[l]->N[d]->v,L->x[d]);
	}
}
/*  */ 

}







/*  updates all the variables in a single chain */
/*  I do this in the order of W and Z first, then theta and */
/*  pi, because it is easier to initialize Pi and Theta (and the Z's) */
/*  than initializing W's and Z's */
void SingleSweep(hyb_chain *C)
{
	
	
	/*  Do all the updates specific for each type of simulation */
	if(C->Lat->TypeOfSim == NEW_HYBRIDS) {
		UpdateWandZandY(C);
		UpdatePi(C);
	}	
	
	
	
	/*  here, this one gets done always */
	UpdateTheta(C);
	
	/*  compute the Kullback-Leibler divergence between the populations */
	KullbLeib(C);
	
	/* do the cline updates if indicated */
	if(gClines->NumLocales) {
		UpdateClineVariables(C, gClines);
	}
	
	
	/* now print the trace if requested  */
	if(C->Dat->PiTraceReport>0 && C->Lat->Pi[0]->NumAved % C->Dat->PiTraceReport == 0) {
		printPi_Trace(C->Lat->Pi[0]->NumAved, C,0); 
	}
	
	
	
}



/*  adds the current realization to the averages and standard errors and this histograms */
/*  of all the variables of interest. */
/*  If DoThetas = 0, then it won't increment the thetas.  This will be useful when  */
/*  we start doing multiple sampling locations */
void IncrementValues(hyb_chain *C, int DoThetas)
{
	int a,l,k,g,i;
	
	/*  increment the thetas */
	if(DoThetas)  {
		CYCLE_alk(C->Dat)
			IncrementDval(C->Lat->Theta[a][l][k]);
		END3CYCLE
	}
	
	/*  increment pi */
	CYCLE_g(C->Lat) 
		IncrementDval(C->Lat->Pi[g]);
	END1CYCLE
	
	/*  now cycle over the individuals and take care of PofZ  */
	/*  in the future I will take care of W, PofW, and the unobserved */
	/*  values for dominant marker types... */
	CYCLE_i(C->Dat)
		CYCLE_g(C->Lat)
			IncrementDval(C->Lat->Ind[i]->PofZ[g]);
			IncrementDval(C->Lat->Ind[i]->UniPriPofZ[g]);
		END1CYCLE
		
		CYCLE_l(C->Dat)  /* a quick loop to do the origins of the particular gene copies */
			for(a=0;a<2;a++) {
				IncrementIval(C->Lat->Ind[i]->W[l][a]);
			}
		END1CYCLE

	END1CYCLE
	
	CYCLE_l(C->Dat)
		IncrementDval(C->Lat->Locus_KB[l]);
	END1CYCLE
	
	IncrementClineVars(C, gClines);
	
}

/*  allocate memory to a hyb_indiv_latent struct, and set the pointers */
/*  for the allelic types to the observed data for codominant locus types. */
/*  It needs the i to set the pointers correctly to the i-th individual */
hyb_indiv *CreateIndiv(hyb_data *D, int i)
{
	int l,w1;
	int w[2],y[2];
	hyb_indiv *temp;
	
	/*  allocate to temp */
	temp = (hyb_indiv *)ECA_MALLOC(sizeof(hyb_indiv));
	
	
	/*  allocate pointers and/or memory to the observed data Y, and to W while we are at it */
	/*  and to PofW as well */
	temp->Y = (ival ***)ECA_MALLOC(D->L * sizeof(ival **));
	temp->W = (ival ***)ECA_MALLOC(D->L * sizeof(ival **));

/*  if(Ind->WsGetSampled==1)  { */
	temp->PofW = (dval ****)ECA_MALLOC(D->L * sizeof(dval ***));
	temp->PofWandYlat = (double *****)ECA_MALLOC(D->L * sizeof(double ****));
/*  } */
/*  else "If there are AFLP's in the data, allocate to temp->PofYlatGivenW" */
/*  this will hold the values for indivs from fixed sources for drawing new Y' */


	/*  set the index of the individual, too */
	temp->idx = i;
	
	CYCLE_l(D)  /*  cycle over the loci */
		temp->Y[l] = (ival **)ECA_MALLOC(2 * sizeof(ival *));
		temp->W[l] = IvalVector(0,1, 0,0,-1);
		
		temp->PofW[l] = (dval ***)ECA_MALLOC(2*sizeof(dval **));
		for(w1=0;w1<2;w1++)  {
			temp->PofW[l][w1] = DvalVector(0,1, 0.0,0.0,-1.0);
		}
		temp->PofWandYlat[l] = NULL;  /*  set it originally to NULL, which is appropriate for  */
								   /*  CODOM markers */
								   
		if(D->LocTypes[l] == AFLP)  {  /*  but if the locus is AFLP, then we must allocate space */
									   /*  to the thing.  This is going to look a little ugly. */
/*  if(Ind->WsGetSampled==1)  {			 */
			/*  For the joint probs of Y and W: */
			temp->PofWandYlat[l] = (double ****)ECA_MALLOC(2 * sizeof(double ***));
			CYCLE_TO_1(w[0])
				temp->PofWandYlat[l][w[0]] = (double ***)ECA_MALLOC(2 * sizeof(double **));
			CYCLE_TO_1(w[1])
				temp->PofWandYlat[l][w[0]][w[1]] = (double **)ECA_MALLOC(2 * sizeof(double *));
			CYCLE_TO_1(y[0])
				temp->PofWandYlat[l][w[0]][w[1]][y[0]] = (double *)ECA_CALLOC(2,sizeof(double));
			END3CYCLE	
/*  } */
/* else  "allocate to a smaller set of PofYsGivenW"   //} */
			/*  For the actual Y value at this locus */
			temp->Y[l] = IvalVector(0,1, 0, 0, 0);					   
			
		}
		
	END1CYCLE
	
	/*  Then cycle over the gtyp freq categories and allocate space for PofZ and LogPofY */
	temp->LogPofY = DvalVector(0,D->Gn->v-1, 0.0, 0.0, -1.0);
	temp->PofZ = DvalVector(0,D->Gn->v-1, 0.0, 1.0, .01);
	temp->UniPriPofZ = DvalVector(0,D->Gn->v-1, 0.0, 1.0, .01);
	
	/*  and allocate space to Z too */
	temp->Z = AllocIval(0,0,0);
	
	/* and point the indiv's flag pointer to the appropriate IndFlags in the Data struct */
	temp->Flag = D->IndFlags[i];
	
	
	/*  And, now for the exciting part, we point all the observed allele pointers to */
	/*  the values in the data struct for loci of type CODOM.   */
	/*  However, for the AFLP loci we don't do that, since these already point to memory */
	/*  that got allocated from them above.    AND!! IF THE LOCUS WAS NOT OBSERVED */
	/*  I.E. WAS MISSING!! THEN WE HAVE TO ALLOCATE MEMORY TO Y.  THEN WE NEED TO WRITE */
	/*  ANOTHER FUNCTION THAT IMPUTES THE VALUES OF Y IF Yobs == -1.  THIS IMPUTATION WILL */
	/*  BE BASED ON THE CURRENT VALUE OF z, I BELIEVE. */
	
	CYCLE_l(D)  /*  cycle over the loci */
		if(D->LocTypes[l] == CODOM) {
			for(w1=0;w1<2;w1++)  {
				temp->Y[l][w1] = D->Yobs[i][l][w1];
			}
		}
	END1CYCLE
	
	
	return(temp);
	
}




/*  Allocates all the necessary memory for the latent data associated with one chain */
/*  it uses the information in the data struct to determine what is necessary */
hyb_chain *CreateLatentChain(hyb_data *D, hyb_prior *P)
{	
	int a,l,i;
	hyb_chain *C;
	
	/*  allocate to C */
	C = (hyb_chain *)ECA_MALLOC(sizeof(hyb_chain));
	
	/*  first thing:  make the data pointer in C point to D */
	C->Dat = D;
	
	/*  then make the prior pointer point to the Priors */
	C->Pri = P;
	
	/*  Now allocate to the latent struct */
	C->Lat = (hyb_latent *)ECA_MALLOC(sizeof(hyb_latent));
	
	/*  make the Gn and G pointers in Lat point to the appropriate places in D */
	C->Lat->Gn = D->Gn;
	C->Lat->G  = D->G;
	
	/*  allocate to Pi array and the s-array */
	C->Lat->Pi = DvalVector(0,C->Lat->Gn->v-1, 0.0, 1.0, .01);
	C->Lat->s = (double *)ECA_CALLOC((size_t)C->Lat->Gn->v,sizeof(double));
	
	/*  by default, make the type of simulation NewHybrids */
	C->Lat->TypeOfSim = NEW_HYBRIDS;
	
	/*  allocate to Theta arrays and r->arrays */
	C->Lat->Theta = (dval ****)ECA_MALLOC(2*sizeof(dval ***));
	C->Lat->r = (double ***)ECA_MALLOC(2*sizeof(double **));
	CYCLE_a
		C->Lat->Theta[a] = (dval ***)ECA_MALLOC(C->Dat->L * sizeof(dval **));
		C->Lat->r[a] = (double **)ECA_MALLOC(C->Dat->L * sizeof(double *));
	END1CYCLE
	
	CYCLE_al(C->Dat)
		C->Lat->Theta[a][l] = DvalVector(0,C->Dat->Kl[l]-1, 0.0, 1.0, .01);
		C->Lat->r[a][l] = (double *)ECA_CALLOC((size_t)C->Dat->Kl[l],sizeof(double));
		
	END2CYCLE
	
	/*  allocate to the array of pointers to individuals */
	C->Lat->Ind = (hyb_indiv **)ECA_MALLOC(C->Dat->M * sizeof(hyb_indiv));
	
	/*  then cycle over the individuals and allocate space to each one of them */
	CYCLE_i(C->Dat)
		C->Lat->Ind[i] = CreateIndiv(D,i);
	END1CYCLE
	
	/*  and allocate space for the Kullback-Leibler divergence of each locus */
	C->Lat->Locus_KB = DvalVector(0,C->Dat->L-1, 0.0, 0.0, 0.0);
	
	return(C);

}


/*  Create Jeffrey's Priors for the Pi and Theta structs */
/*  This returns a pointer to the priors and takes as */
/*  input the Data, D */
hyb_prior * CreatePriors(hyb_data *D, enum prior_type PiPriorType, enum prior_type ThetaPriorType)
{
	int l,a;
	hyb_prior *P;
	
	/*  first allocate the necessary memory */
	/*  To Lambda (prior parameters for Theta) */
	P = (hyb_prior *)ECA_MALLOC(sizeof(hyb_prior));
	P->Lambda  = (dval ****)ECA_MALLOC(2*sizeof(dval ***) );
	CYCLE_a
		P->Lambda[a] = (dval ***)ECA_MALLOC(D->L * sizeof(dval**));
		CYCLE_l(D)
			P->Lambda[a][l] = DvalVector(0,D->Kl[l]-1, 0.0,0.0,-1.0); 
		END1CYCLE
	END1CYCLE
	/*  to Zeta (prior parameters for Pi) */
	P->Zeta = DvalVector(0,D->Gn->v-1, 0.0,0.0,-1.0);
	
	/*  then set the values */
	SetPiPriors(D,P,PiPriorType);
	SetThetaPriors(D,P,ThetaPriorType);
	
	/*  may as well allocate and set values for the JABES related variable */
	P->XiPrior = DvalVector(0,1,0.0,1.0,.02);
	/*  and set them to Jeffreys right now too: */
	P->XiPrior[0]->v = .5;
	P->XiPrior[1]->v = .5;
	
	return(P);
}



/*  puts  Prior values in for the mixing proportions */
/*  assumes that memory has already been allocated to Zeta array */
void SetPiPriors(hyb_data *D, hyb_prior *P, enum prior_type PiPriorType)
{
	int g;
	
	CYCLE_g(D)
		if(PiPriorType == JEFFREYS)
			P->Zeta[g]->v = 1.0/(double)(D->Gn->v);
		else if(PiPriorType == UNIFORM)
			P->Zeta[g]->v = 1.0;
		else  /*  the default is Jeffreys */
			P->Zeta[g]->v = 1.0/(double)(D->Gn->v);
	END1CYCLE
	
/*	if(PiPriorType == JEFFREYS)
		P->PiPriType = JEFFREYS;
	else if(PiPriorType == UNIFORM)
		P->PiPriType = UNIFORM;
	else 
		P->PiPriType = OTHER;
*/
	P->PiPriType = PiPriorType;
}

/*  puts  Prior values in for the mixing proportions */
/*  assumes that memory has already been allocated to Lambda arrays */
void SetThetaPriors(hyb_data *D, hyb_prior *P, enum prior_type ThetaPriorType)
{
	int a,l,k;
	
	CYCLE_alk(D)
		if(ThetaPriorType == JEFFREYS)
			P->Lambda[a][l][k]->v = 1.0/(double)(D->Kl[l]);
		else if(ThetaPriorType == UNIFORM)
			P->Lambda[a][l][k]->v = 1.0;
		else  /*  default is Jeffreys prior */
			P->Lambda[a][l][k]->v = 1.0/(double)(D->Kl[l]);
	END3CYCLE
	
	if(ThetaPriorType == JEFFREYS)
		P->ThetaPriType = JEFFREYS;
	else if(ThetaPriorType == UNIFORM)
		P->ThetaPriType = UNIFORM;
	else 
		P->ThetaPriType = USER_DIRCH;
}


/*  Initialize Gn, Theta, Pi, (and Z) to "overdispersed" starting values */
/*  by simulating Theta and Pi from Dirichlet Dsns given by their priors */
/*  and then simulating Z's from the resulting Pi. */
/*  THIS USE_FIXED_ZS_4_PI_INITSHOULD ALSO INITIALIZE Ind->Y TO RANDOM ALLELIC VALUES FOR LOCI THAT */
/*  ARE MISSING DATA. */
void InitializeChain(hyb_chain *C)
{
	int a,l,i,g,k;
	double tempTheta[MAX_ALLELES], tempPrior[MAX_ALLELES];
	double tempPi[MAX_GN], tempPiPrior[MAX_GN];
	
	
	/*  YO!  FIGURE OUT WHAT VALUE YOU WILL USE TO SET GN */
	/*  TO SOMETHING HERE. */
	
	
	/*  cycle over species and loci, initializing Theta */
	CYCLE_al(C->Dat) 
		/*  store the priors in the temporary array */
		CYCLE_k(C->Dat)
			tempPrior[k] = C->Pri->Lambda[a][l][k]->v;
			if(C->Dat->StartingDispersionTheta==USE_PRIORIZED_ALLELES)  {
				tempPrior[k] += C->Dat->PriorizedAlleleCounts[a][l][k];
			}
		END1CYCLE
		/*  simulate the new value */
		DirichletRV(tempPrior,C->Dat->Kl[l],tempTheta);
		
		/*  copy it into the theta struct */
		CYCLE_k(C->Dat)
			C->Lat->Theta[a][l][k]->v = tempTheta[k];
		END1CYCLE
	END2CYCLE
	
	
	
	/*  then put the Pi Priors into temp storage and draw the initial Pi */
	CYCLE_g(C->Lat)
		tempPiPrior[g] = C->Pri->Zeta[g]->v;
	END1CYCLE
	
	if(C->Dat->StartingDispersionPi==USE_FIXED_ZS_4_PI_INIT) {
		CYCLE_i(C->Dat)
			if(C->Lat->Ind[i]->Flag[Z_FIXED])  {
				if(C->Lat->Ind[i]->Flag[DONT_CONTRIBUTE_TO_PI]==0) {
					tempPiPrior[C->Lat->Ind[i]->Flag[FIXED_Z]] += 1.0;
				}
			}
		END1CYCLE
	}
	DirichletRV(tempPiPrior,C->Lat->Gn->v,tempPi);
			
	/*  then copy those Pi values into the appropriate places */
	CYCLE_g(C->Lat)
		C->Lat->Pi[g]->v = tempPi[g];
	END1CYCLE
	
	/*  then simulate Z for each individual by drawing directly from Pi.
		Unless Z is fixed, in which case you just fix Z as it should be */
	CYCLE_i(C->Dat)
		if(C->Lat->Ind[i]->Flag[Z_FIXED]==1) 
			C->Lat->Ind[i]->Z->v = C->Lat->Ind[i]->Flag[FIXED_Z];
		else
			C->Lat->Ind[i]->Z->v = IntFromProbsRV(tempPi,0,C->Lat->Gn->v);
		
		/* the two following procedures suggested in the comments are not implemented,
		because I don't think it would really make much of a difference  */
			/*  then set the W's randomly given Z for each individual. */
			
			/*  Also, set the W values for the known source individuals who are */
			/*  not getting prior-ized.  Even for the AFLP ones. */
	END1CYCLE
	
	
	RandoStartClineVars(C, gClines);
	
	/* and here we can print headers for trace reports that we will be doing */
	if(C->Dat->PiTraceReport>0) {
		printPi_Trace(0, C, 1);
	}
}





/*  determines the value of r---the number of each type of allele */
/*  currently allocated to one or the other of the populations.  This */
/*  is currently written with codominant loci in mind */
/*  THIS WILL HAVE TO BE CHANGED FOR OTHER TYPES OF LOCUS */
void Count_r(hyb_chain *C)
{
	int i,l,k,a,c;  /*   subscripts for individuals (i),loci (l), alleles (k),  */
				  /*   species (a) or gene-copies (c) */
	
	/*  first, initialize r-values to 0's */
	CYCLE_alk(C->Dat)
		C->Lat->r[a][l][k] = 0.0;
	END3CYCLE
	
	/*  then cycle over the individuals and the gene copies within them */
	/*  and count up the alleles which are allocated according to W */
	CYCLE_i(C->Dat)
		CYCLE_l(C->Dat)
			
			/*  only count up the non-missing alleles.  Test here before looking at each allele */
			/*  (if one is missing, the other is too!) 
				WITH FIXED_Z OPTIONS:  don't add in the alleles from CODOM loci possessed by FIXED_IN_PURE
				individuals, because those have already been incorporated into the allele frequency priors.  */
			if(YnotMissing(C,i,l) && !(C->Lat->Ind[i]->Flag[Z_FIXED_IN_PURE] && C->Dat->LocTypes[l] == CODOM)    \
			        && !(C->Lat->Ind[i]->Flag[DONT_CONTRIBUTE_TO_THETA]) ) {
				CYCLE_c
					C->Lat->r[ C->Lat->Ind[i]->W[l][c]->v ]
							 [ l ]
							 [ C->Lat->Ind[i]->Y[l][c]->v ]  += 1.0;
				END1CYCLE
			}
	END2CYCLE
}



/*  Determines the value of s---the number of individuals currently allocated to */
/*  each of the different hybrid classes */
void Count_s(hyb_chain *C)
{
	int g,i;
	
	/*  first, initialize s values to 0.0's */
	CYCLE_g(C->Lat)
		C->Lat->s[g] = 0.0;
	END1CYCLE
	
	/*  then cycle over the individuals and add 1.0's to s according to their Z's */
	CYCLE_i(C->Dat)
		if(C->Lat->Ind[i]->Flag[DONT_CONTRIBUTE_TO_PI] != 1)
			C->Lat->s[  C->Lat->Ind[i]->Z->v ] += 1.0;
	END1CYCLE

}


/*  a single function that updates all the theta values */
void UpdateTheta(hyb_chain *C)
{
	int a,l,k;
	double temp[MAX_ALLELES];  /*  a temporary array of doubles */
	
	/*  first count up the r's */
	Count_r(C);
	
	/*  then add the priors to those and the "priorized" alleles */
	CYCLE_alk(C->Dat)
		C->Lat->r[a][l][k] += C->Pri->Lambda[a][l][k]->v + C->Dat->PriorizedAlleleCounts[a][l][k];
	END3CYCLE
	
	/*  Then cycle over pops and loci, simulating new values into the temporary vector */
	/*  and then copying them into the Theta structs */
	CYCLE_al(C->Dat)
		/*  simulate the new value */
		DirichletRV(C->Lat->r[a][l],C->Dat->Kl[l],temp);
		
		/*  copy it into the theta struct */
		CYCLE_k(C->Dat)
			C->Lat->Theta[a][l][k]->v = temp[k];
		END1CYCLE
	END2CYCLE
}


/*  a single function that updates all the pi values */
void UpdatePi(hyb_chain *C)
{
	int g;
	double temp[MAX_GN];  /*  a temporary array of doubles */
	
	if(gPiFixed==1) {
		CYCLE_g(C->Lat)
			C->Lat->Pi[g]->v = gPiFixedValues[g];
		END1CYCLE
		return;
	}
	else {
		/*  first count up the r's */
		Count_s(C);
		
		/*  then add the priors to those */
		CYCLE_g(C->Lat)
			C->Lat->s[g] += C->Pri->Zeta[g]->v;
		END1CYCLE
		
		/*  then simulate a new value  */
		DirichletRV(C->Lat->s,C->Lat->Gn->v,temp);
		
		/*  and copy it into the Pi struct */
		CYCLE_g(C->Lat)
			C->Lat->Pi[g]->v = temp[g];
		END1CYCLE
	}
	
}


/*  for a single individual, this computes and stores the Log of the probabilities */
/*  of the genotypes of the individual for all the possible Gn gtyp freq categories */
/*  Along the way, it stores the Probability of the W's for that individual given */
/*  its current gtyp freq category allocation (i.e. given its current Z). */
/*  We pass in C so that we know the limits to loop over (i.e. the number of loci, etc). */
double StoreSingleLogPofY(hyb_indiv *Ind, hyb_chain *C)
{	
	int g,l;

	/*  cycle over gtyp freq classes */
/*  if(Ind->Flag->LogPofY==1 || C->Flag->DOM_Present == 1)  */
	CYCLE_g(C->Lat)
	
/*  if(Ind->Flag->LogPofY==1 || Ind->ZisFixed) */
		/*  initialize the LogPofY[g] to 0.0 to accumulate a sum over loci */
		Ind->LogPofY[g]->v = 0.0;
		
		/*  cycle over loci */
		CYCLE_l(C->Dat)
			Ind->LogPofY[g]->v += log(PofY_at_a_locus(Ind,C,l,g));
		END1CYCLE  /*  l */
		
	END1CYCLE /*  g */
/*  } */
	return(Ind->LogPofY[Ind->Z->v]->v);  /*  return the prob of its gtyp given the class it is  */
										 /*  allocated to. */
}


/*  computes the probability of Y at a locus l conditional on the individual */
/*  belonging to gtyp frequency class g.  Basically this implements */
/*  Equation 2 in the paper, summing, but summing over the four possible W's */
/*  C gives us a way to access the array G of expected proportions of the different genotypes */
/*  for the g-th genotype frequency class. */
/*  Crucially, this function also stores the full conditionals for W in PofW, based on the current */
/*  z of the individual. */
/*  For an AFLP locus, it stores the appropriate probabilities of the pairs (W,Y) because */
/*  the Y's are latent in the AFLP case.    */
double PofY_at_a_locus(hyb_indiv *Ind,
		hyb_chain *C, int l, int g)
{
	int i;
	double temp;	/*  initialize to accumulate a product */
	double normo=0.0;  /*  initialize to accumulate a sum */
	int w[2], y[2];
	
	
	/*  	DEALING WITH MISSING DATA: */
	/*  for now, we are not imputing, so if the data are missing, we */
	/*  just return 1.0, whether or not the locus is Codominant or AFLP */
	if(!YnotMissing(C,Ind->idx,l))
		return(1.0);
		
	
	
	
	
	/*  CODOM ROUTINE */
	if(C->Dat->LocTypes[l] == CODOM)  {
		DO_W_CYCLE(w)   /*  cycle over the different possible W's to marginalize them out */
			temp = 1.0;  /*  initialize to accumulate a product for each w[0],w[1] pair. */
			
			/*  get the allele frequency part in there first */
			for(i=0;i<2;i++)  {
				temp *= C->Lat->Theta[ w[i] ][l][ Ind->Y[l][i]->v ]->v;
			}
			
			/*  multiply the G part in there.  This depends on the gtyp freq class */
			/*  g that we are interested in */
			temp *= C->Lat->G[g][ w[0] ][ w[1] ];
			
			normo += temp;
					
			/*  now deal with storing the PofW part for the locus */
			if(g==Ind->Z->v)  {
				Ind->PofW[l][ w[0] ][ w[1] ]->v = temp;
			}
		END2CYCLE
		
		/*  now, if we just did the g for which the individual is currently allocated */
		/*  then we need to cycle over the w's one more time to normalize the PofW */
		if(g==Ind->Z->v)  {
			DO_W_CYCLE(w)
				Ind->PofW[l][ w[0] ][ w[1] ]->v /= normo;
			END2CYCLE
		}
	} 
	
	
	/*  AFLP ROUTINE:  For all g except the current one, we compute the probability of Y  */
	/*  given it is a "-" */
	/*  and then we either use that, if it really is a "-" or subtract it from 1.0 if */
	/*  it really is a "+".  For g that is the current one, we sum over all the possibilities */
	/*  so we can store them and use them to draw the W and the latent Y later.  Note that */
	/*  the "+" allele will always be allele 1 and the "-" allele will always be allele 0.  */
	if(C->Dat->LocTypes[l] == AFLP)  {

		if(g != Ind->Z->v)  {
			DO_W_CYCLE(w)   /*  cycle over the different possible W's to marginalize them out */
							/*  conditional on the whole individual being "-" for the locus */
				temp = 1.0;  /*  initialize to accumulate a product for each w[0],w[1] pair. */
				
				/*  get the allele frequency part in there first */
				for(i=0;i<2;i++)  {
					temp *= C->Lat->Theta[ w[i] ][l][ 0 ]->v;   /*  "0" subscripts the "-" allele */
				}
				
				/*  multiply the G part in there.  This depends on the gtyp freq class */
				/*  g that we are interested in */
				temp *= C->Lat->G[g][ w[0] ][ w[1] ];
				
				normo += temp;
			END2CYCLE
			
			/*  now, if the individual is really a "-", then normo is the probability that */
			/*  we want.  But if not, then we want to reset normo to 1-normo; */
			/*  Note that Yobs[i][l][0] (NOT Yobs[i][l][1]) is where the 0 or 1, */
			/*  signifying a "-" or a "+" resides. */
			if(C->Dat->Yobs[Ind->idx][l][0]->v == 1) 
				normo = 1.0 - normo; 
		}
		
		
		/*  if g is the the class this individual is currently allocated to, then we have to cycle */
		/*  over both W and possible latent y's, and store all that stuff. */
		else {
			
			DO_W_CYCLE(w)
				DO_W_CYCLE(y)
					
					/*  compute the probs of the w and y pairs.  if it is going to */
					/*  be zero, evaluate that first, to short-circuit multiplication by the */
					/*  allele freqs.  Basically if either of y[0] or y[1] are 1 when the  */
					/*  observed state is "-"  (coded as zero) then the prob is zero.  The same */
					/*  goes for y[0] and y[1] both being 0 when the observed state is "+" */
					if(   (  (y[0] + y[1] > 0) && (C->Dat->Yobs[Ind->idx][l][0]->v == 0 ) )  ||
						  (  (y[0] + y[1] == 0) && (C->Dat->Yobs[Ind->idx][l][0]->v == 1 ) )	)
					{
						temp = 0.0;
					}
					else {
						/*  first put the G part in there */
						temp = C->Lat->G[g][ w[0] ][ w[1] ];
						
						for(i=0;i<2;i++)  {
							temp *= C->Lat->Theta[ w[i] ][l][ y[i] ]->v;
						}
					}
					
					normo += temp;
					Ind->PofWandYlat[l][w[0]][w[1]][y[0]][y[1]] = temp;
				END2CYCLE
			END2CYCLE
			
			/*  finally, cycle one last time to normalize the things */
			DO_W_CYCLE(w)
				DO_W_CYCLE(y)
					Ind->PofWandYlat[l][w[0]][w[1]][y[0]][y[1]] /= normo;
				END2CYCLE
			END2CYCLE
			
		}
	}
	
	
	/*  that's it then! */
	return(normo);
}



/*  updates the W's and Z's for all the individuals, and the Latent Y's */
/*  for the AFLP loci.   */
/* here define WAY_TOO_SMALL to the the log of 10 to the -200 */
#define WAY_TOO_SMALL -460.517
void UpdateWandZandY(hyb_chain *C) 
{
	int i,g;
	double normo, uniprinormo, log_extract=-99999999999999999999.0, tmp_prob;
	double TempCompleteDataLogLike = 0.0;  /*  prepare to accumulate a sum */
	
	CYCLE_i(C->Dat)
		/*  compute the genotype probs (and the full conditionals for the w's) */
		/*  for each individual */
		TempCompleteDataLogLike += StoreSingleLogPofY(C->Lat->Ind[i], C);
		
		/*  then we must apply Bayes law to use PofY and Pi to compute  */
		/*  PofZ.  In the future I may want to do some fancy stuff to prevent */
		/*  underflow here (i.e. add in and later pull out a constant of log something before */
		/*  something before exponentiating LogPofY to become a probability).  */
		
		normo = 0.0;  /*  initialize to accumulate a sum of probabilities */
		uniprinormo = 0.0; 

/*printf("Ind %d : ",i);*/
    
    /*  cycle once over the genotype frequency classes to find out which is the max log-prob */
    CYCLE_g(C->Lat)
      if(C->Lat->Ind[i]->LogPofY[g]->v > log_extract) log_extract = C->Lat->Ind[i]->LogPofY[g]->v;
    END1CYCLE
  

		/*  cycle over the genotype frequency classes */
		CYCLE_g(C->Lat)
      
    if(C->Lat->Ind[i]->LogPofY[g]->v - log_extract < WAY_TOO_SMALL)
      tmp_prob = 0.0;
    else 
      tmp_prob = exp(C->Lat->Ind[i]->LogPofY[g]->v - log_extract);
    

    
			/*  first store them as unnormalized probabilities */
			C->Lat->Ind[i]->PofZ[g]->v = tmp_prob * C->Lat->Pi[g]->v;


/*printf(" %.4f",C->Lat->Ind[i]->LogPofY[g]->v);*/

			/* store the scaled likelihoods too */				 
			C->Lat->Ind[i]->UniPriPofZ[g]->v = tmp_prob;
										 
			/*  and add that to the normalization factor */
			normo += C->Lat->Ind[i]->PofZ[g]->v;
			uniprinormo += C->Lat->Ind[i]->UniPriPofZ[g]->v;
			
		END1CYCLE
/*printf("\n");*/
#undef WAY_TOO_SMALL
  
  
		/*  now cycle again over g and normalize PofZ */
		CYCLE_g(C->Lat)
			/*  divide each one by normo */
			C->Lat->Ind[i]->PofZ[g]->v /= normo;
			C->Lat->Ind[i]->UniPriPofZ[g]->v /= uniprinormo;
		END1CYCLE


/* Draw the Z first, and then recompute the probabilities for the W's (this way it
will store the W's for the new Z, and use those to choose the new W's.  I used to do this
in a silly order---choose the w's then choose the z.  This caused the program to occasionally
get caught up in fast and furious label switching. */
if(C->Lat->Ind[i]->Flag[Z_FIXED] == 0)  {  /* only do this for non-fixed-Z individuals */
	C->Lat->Ind[i]->Z->v = IntegerFromProbsRV(C->Lat->Gn->v,C->Lat->Ind[i]->PofZ);
	StoreSingleLogPofY(C->Lat->Ind[i], C); /* re-compute the new W probs for this new Z */		
}
		
		/*  now we are set to draw new values for the W's and the Z's for each individual */
		/*  draw the W's first */
		DrawNewW_and_LatentY_indiv(C->Lat->Ind[i],C);
		
		
// THIS BIT WAS COMMENTED OUT, IT IS SUPERSEDED NOW BY THE STUFF 7 to 14 lines ABOVE		
		/*  then draw the Z for the individual (you CANNOT do this before */
		/*  drawing the new W, without screwing things up!! */
//		if(C->Lat->Ind[i]->Flag[Z_FIXED] == 0)  {  /* only do this for non-fixed-Z individuals */
//			C->Lat->Ind[i]->Z->v = IntegerFromProbsRV(C->Lat->Gn->v,C->Lat->Ind[i]->PofZ);		
//		}
		
	END1CYCLE
	
	/*  now, at the end of all that, we update the TempCompleteDataLogLike Trace */
	C->Lat->TheCompleteDataLogLike = TempCompleteDataLogLike;
}



/*  returns a discrete random integer between 0 and Num-1  */
/*  based on probabilities stored in an array of pointers to dvals. */
/*  This assumes that the probabilities are subscripted 0,...,Num-1 */
/*  Does not check to make sure that the probabilities */
/*  sum to one.  Requires function ranf() from ranlib or elsewhere */
/*  to generate a Uniform double on (0,1). */
/*  This should go into the MathStatRand library, eventually  */
int IntegerFromProbsRV(int Num,dval **P)
{
	int i;
	double cum=0.0;
	double rando = (double)ranf();
	
	for(i=0;i<Num;i++)  {
		cum += P[i]->v;
		if(cum > rando)
			break;
	}
	if(i==Num && cum<.9999999999)  {
  		fprintf(stderr,"\n\nCumulative Prob less than 1.0 in IntegerFromProbsRV\n\n");
		fprintf(stderr,"Exiting to system...\n\n");
		exit(1);
	}
	
	/* now, sometimes if the P->v's are not well defined or are NaN's things get
	   screwed up and we end up at this point with i==Num.  (This is a horrible thing).
	  	I have seen it happen on the first sweep, as if things haven't gotten initialized
	  	just right, but then, if the value gets restored to the right on it does OK.  So,
	  	I am going to fix that here by just returning a uniform int and barking a warning to 
	  	stderr each time it happens. */
	
	if(i>=Num) {
		fprintf(stderr,"WARNING: i=%d >= Num=%d in IntegerFromProbsRV().  Probably only a problem while initializing the chain. ",i,Num);
		i = UniformRV(0,Num-1);
		
		fprintf(stderr," i reset to %d\n",i);
	}
	
	return(i);
}




/*  cycle over the all the loci in an individual and  */
/*  update the W values at each locus.  This requires that */
/*  PofW be computed previously (usually by a call to  */
/*  JointPofWandY_at_a_locus) and that the Z for the individual */
/*  has not been changed since then.   */
void DrawNewW_and_LatentY_indiv(hyb_indiv *Ind, hyb_chain *C)
{
	int l,w[2],y[2];
	double cum,rando;
	
	CYCLE_l(C->Dat)
	
		/*  ONLY DRAW A NEW W and LATENT Y if THE LOCUS IS NOT MISSING!!! */
		if(YnotMissing(C,Ind->idx,l)) {
			cum = 0.0;  /*  set to zero to accumulate a sum */
			rando = (double)ranf();
	/* if(Ind->SampleWs==1)  {		 */
			/*  CODOM ROUTINE: */
			if(C->Dat->LocTypes[l] == CODOM)  {
				DO_W_CYCLE(w) 
					cum += Ind->PofW[l][ w[0] ][ w[1] ]->v;
					if(cum>rando) {
						/*  set the two values of Ind->W accordingly */
						Ind->W[l][0]->v = w[0];
						Ind->W[l][1]->v = w[1];
						/*  then set the little w's each to 2 so that we exit the two for loops */
						w[0] = 2;
						w[1] = 2;
					}
				END2CYCLE
			}
			
			/*  AFLP ROUTINE */
			if(C->Dat->LocTypes[l] == AFLP)  {
				DO_W_CYCLE(w)
		 			DO_W_CYCLE(y)
						cum += Ind->PofWandYlat[l][ w[0] ][ w[1] ][ y[0] ][ y[1] ];
						if(cum>rando) {
							/*  set the two values of Ind->W and of Ind->Y accordingly */
							Ind->W[l][0]->v = w[0];
							Ind->W[l][1]->v = w[1];
							Ind->Y[l][0]->v = y[0];
							Ind->Y[l][1]->v = y[1];
							/*  then set the little w's and y's all to 2 so that we exit the two for loops */
							w[0] = 2;
							w[1] = 2;
							y[0] = 2;
							y[1] = 2;
						}
					END2CYCLE
				END2CYCLE
			}
	/*  }   // closes if(Ind->SampleWs==1) */

	/*  else if(C->Dat->LocTypes[l] == AFLP)  "Sample The Ys but not Ws at This AFLP Locus"  } */
		}
	
	END1CYCLE
}



/*  Writes out the histograms to a series of files */
void OutputHistograms(hyb_chain *C)
{
	int l,k,a,g;
	FILE *out;
	char FILESTRING[10000];
	
	sprintf(FILESTRING,"%saa-Theta.hist",gPWD);
	
	/*  deal first with the thetas: */
	fprint_HistTopRow(C->Lat->Theta[0][0][0]->Hist, FILESTRING);
	out = fopen(FILESTRING,"a");
	
	CYCLE_lka(C->Dat)
		fprintf(out,"\nLoc%d_Alle%d_Pop%d  ",l,k,a);
		fprint_HistLineProp(out,C->Lat->Theta[a][l][k]->Hist);
	END3CYCLE
	fclose(out);
	
	/*  deal next with the Pi's */
	sprintf(FILESTRING,"%saa-Pi.hist",gPWD);
	fprint_HistTopRow(C->Lat->Theta[0][0][0]->Hist, FILESTRING);
	out = fopen(FILESTRING,"a");
	
	CYCLE_g(C->Lat) 
		fprintf(out,"\n%.5f/%.5f/%.5f  ",C->Lat->G[g][0][0],
									  C->Lat->G[g][0][1] + C->Lat->G[g][1][0],
									  C->Lat->G[g][1][1]);
		fprint_HistLineProp(out,C->Lat->Pi[g]->Hist);
	END1CYCLE
	fclose(out);
}






/*  converts a locus_type enum to a string */
const char *LocType2str(enum locus_type loc)
{

	switch (loc) {
	case CODOM: return "CODOMINANT";
	case AFLP: return "AFLP";
	case UNSET: return "UNSET";
	case HAS_NULL: return "HAS_NULL";
	default: return "!!Error-Unknown Locus Type!!";
	}
}





/*  Reset all the averages and things for all of the variables in a chain */
/*  Currently it only deals with the few variables I have been keeping track of */
void ResetAllAveragesEtc(hyb_chain *C)
{
	int a,l,k,g,i;
	/*  first do all the Thetas */
	CYCLE_alk(C->Dat)
		InitDvalSummaryToZero(C->Lat->Theta[a][l][k]);
	END3CYCLE
	
	/*  then do the Pi */
	CYCLE_g(C->Dat) 
		InitDvalSummaryToZero(C->Lat->Pi[g]);
	END1CYCLE

	/*  now cycle over the individuals and take care of PofZ  */
	/*  in the future I will take care of W, PofW, and the unobserved */
	/*  values for dominant marker types... */
	CYCLE_i(C->Dat)
		CYCLE_g(C->Lat)
			InitDvalSummaryToZero(C->Lat->Ind[i]->PofZ[g]);
			InitDvalSummaryToZero(C->Lat->Ind[i]->UniPriPofZ[g]);
		END1CYCLE
		
		CYCLE_l(C->Dat)  /* a quick loop to do the origins of the particular gene copies */
			for(a=0;a<2;a++) {
				InitIvalSummaryToZero(C->Lat->Ind[i]->W[l][a]);
			}
		END1CYCLE

		
	END1CYCLE
	
	
	
	/*  reset the averages of the kullback-leibler stuff too */
	CYCLE_l(C->Dat)
		InitDvalSummaryToZero(C->Lat->Locus_KB[l]);
	END1CYCLE
	
	ResetClineVarAverages(C, gClines);

}


/*  returns 1 if locus l in individual i is not missing */
/*  and  0 if it is missing */
int YnotMissing(hyb_chain *C, int i, int l)
{
	int temp = 1;
	
	if(C->Dat->Yobs[i][l][0]->v == -1 ||  C->Dat->Yobs[i][l][1]->v == -1)  {
		temp = 0;
	}
	
	return(temp);
}



/*  computes the current value of the Kullback-Leibler divergence at a locus */
double KullbLeibAtLocus(hyb_chain *C, int l)
{
	int a,k;
	double temp = 0.0;  /*  initialize to accumulate a sum at this locus */
	
		
	CYCLE_a
		CYCLE_k(C->Dat)
			temp += C->Lat->Theta[a][l][k]->v * log(C->Lat->Theta[a][l][k]->v/
				C->Lat->Theta[(a+1)%2][l][k]->v);  /*  note the little (a+1)%2 trick */
												   /*  to make the denominator the freq */
												   /*  in the other population */
		END1CYCLE
	END1CYCLE
	
	
	/*  now return temp */
	return(temp);
	
	
}


/*  computes all the current Kullback-Leibler divergences at all the loci */
void KullbLeib(hyb_chain *C)
{
	int l;
	
	CYCLE_l(C->Dat)
		C->Lat->Locus_KB[l]->v = KullbLeibAtLocus(C,l);
	END1CYCLE
} 







/*  a silly, expedient function for the scottish cats data set:  starting from */
/*  individual N (numbered N, not its index!!), till the end of the data set, use those individuals to form the */
/*  prior for allele frequencies (along with whatever existing prior had been there), */
/*   and then set the number of individuals to N, so that */
/*  those final ones get dropped from the dataset. */
void UseLastFromNForPriorForZero(hyb_chain *C, int N)
{
	int i,l,c;
	
	for(i=N-1; i<C->Dat->M ; i++)  {
		CYCLE_l(C->Dat)
			if( YnotMissing(C,i,l) )  {
				for(c=0;c<2;c++)  {
					C->Pri->Lambda[0][l][  C->Dat->Yobs[i][l][c]->v ]->v += 1.0;
				}
			}
		END1CYCLE
	}
	
	/*  then, at the end of this we set the number of individuals to N - 1 */
	C->Dat->M = N-1;
		
}


/*  This function counts up all the alleles in the individuals that have fixed
	Z and adds them to the PriorizedAlleleCounts (it also does some memory allocation.
	This is done for CODOM loci.  AFLP loci get PriorizedAlleleCounts of 0 because 
	we don't know if individuals carry one or two copies of the + gene.  */
void PriorizeAllelesFromFixedZ(hyb_data *D)
{
int i,l,c;  /*   subscripts for individuals (i),loci (l), alleles (k),  */
				  /*   species (a) or gene-copies (c) */
	
	/* allocate to PriorizedAlleleCounts using calloc.  This initializes the elements to zero */
	D->PriorizedAlleleCounts = (double ***)ECA_CALLOC((size_t)2, sizeof(double **));
	CYCLE_c
		D->PriorizedAlleleCounts[c] = (double **)ECA_CALLOC((size_t)D->L, sizeof(double *));
		CYCLE_l(D)
			D->PriorizedAlleleCounts[c][l] = (double *)ECA_CALLOC((size_t)D->Kl[l], sizeof(double));
		END1CYCLE
	END1CYCLE
	
	/*  then cycle over the individuals and the gene copies within them */
	/*  and count up the alleles which are allocated according to FIXED_Z */
	CYCLE_i(D)
		if(D->IndFlags[i][Z_FIXED_IN_PURE] == 1) {  /* only do the ones known to come from one of the  
															 two pure gtyp freq cats */
			CYCLE_l(D)
				
				/*  only count up the non-missing alleles.  Test here before looking at each allele */
				/*  (if one is missing, the other is too!).  Also, add counts in only for CODOM loci! */
				if(  !(D->Yobs[i][l][0]->v == -1 ||  D->Yobs[i][l][1]->v == -1) 
						&& D->LocTypes[l] == CODOM ) {
					CYCLE_c
						D->PriorizedAlleleCounts[ D->IndFlags[i][PURE_POOL]]
								 [ l ]
								 [ D->Yobs[i][l][c]->v ]  += 1.0;
					END1CYCLE
				}
			END1CYCLE
		}
	END1CYCLE
}


/* this function asks the user for a file name for the priors and then goes
and adds those to PriorizedAlleleCounts.  This must be called after 
PriorizeAllelesFromFixedZ().  If File is NULL it prompts for a file name
otherwise it uses the one that is there.
*/
void AddPriorizeAllelesFromFile(hyb_data *D, char *InFile)
{
	int temp,lineno=0;
	double temp_d, boost0, boost1;
	char File[300],LocName[100],AlleName[100];
	FILE *in,*out=fopen("aa-ProcessedPriors.txt","w");
	char Str[300];
	int MatchesLocus=0, MatchesAllele=0, loc,alle, get_it=1;
	
	
	if(InFile==NULL) {
	
		/*  now read in the gtyp freq categories */
		printf("\n\n-----------------\n");
		printf("\nEnter the name of a file holding prior allele frequency information");
		printf("\nor just type 0 to not include any prior information\n->");
		
		/*  read the GTYP FILE name or a 0 */
		temp = erdGetNext(File,&temp_d,stdin); 

		
		if(temp == 1 && temp_d == 0.0) {  /*  if input was an int equal to zero just get out*/
			return;
		}
		else {
			in = erdOpenFileOrRetry(File, "r");
		}
	}
	else {
		if( (in=fopen(InFile,"r"))==NULL) {
			fprintf(stderr,"Error! Cannot open file %s to compile information about allele frequency priors. Exiting...\n",InFile);
		}
	}
	
	/* now cycle through the lines of in and match up those that start with "Locus" and after
	  each of those, add the priors in */
	while(!feof(in)) {
		if(get_it==0) 
			get_it=1;
		else {
			fgets(Str,300,in);
			lineno++;
		}
		MatchesLocus = sscanf(Str,"Locus %s",LocName);
		if(MatchesLocus>0)  {
			loc = ReturnNameIndex(D->LocNames,LocName,D->L);
			if(loc>=0)  {
				printf("\nIncorporating Prior for locus %s",LocName);
				fprintf(out,"\nIncorporating Prior for locus %s",LocName);
				/* now we process the alleles---every line that has an int and two doubles...
				   until we hit Locus again or we hit the end of the file  */
				while(!feof(in))  {
					fgets(Str,300,in);
					lineno++;
					if((MatchesAllele = sscanf(Str," %s %lf %lf",AlleName,&boost0,&boost1 ))==3) { /* if it has a line to change things to */
						alle = ReturnNameIndex(D->AlleleNames[loc],AlleName,D->Kl[loc]);
						if(alle>=0) {
							printf("\n\tAllele %s: Add to 0: %f\t\tAdd to 1: %f",AlleName,boost0,boost1);
							fprintf(out,"\n\tAllele %s: Add to 0: %f\t\tAdd to 1: %f",AlleName,boost0,boost1);
							D->PriorizedAlleleCounts[0][loc][alle] += boost0;
							D->PriorizedAlleleCounts[1][loc][alle] += boost1;
						}
						else if(alle==-1) {
							printf("\n\nError!  Allele %s not found at locus %s. Line %d of file \"%s\"\n\nexiting to system...\n\n",AlleName,LocName,lineno,File);
							exit(1);
						}
						else if(alle==-2) {
							printf("\n\nError!  More than two alleles named %s at locus %s. Line %d of file \"%s\"\n\nexiting to system...\n\n",AlleName,LocName,lineno,File);
							exit(1);
						}
						
					} /* closes if(sscanf==3) */
					else if(sscanf(Str,"Locus %s",LocName))  {
						get_it=0;
						break;  /* get out of the current while loop */
					}
				}  /* closes the while loop */
			}	
			else if(loc==-1) {
				printf("\n\nError!  Locus %s not found in the data. Line %d of file \"%s\"\n\nexiting to system...\n\n",LocName,lineno,File);
				exit(1);
			}
			else if(loc==-2) {
				printf("\n\nError!  More than one locus named %s. Line %d of file \"%s\"\n\nexiting to system...\n\n",LocName,lineno,File);
				exit(1);
			}
		}	
	}
	fclose(out);
	fclose(in);
	printf("\n\nTake a moment to review the updates to the priors listed above.");
	printf("\nThe above information is also stored in the file \"aa-ProcessedPriors.txt\"\n");
	
}

/* given a name in Str, this returns the index of a solitary occurrence 
in the first n elements of
Str in array.  If there is more than one occurrence, it returns -2, if there
are no occurrences, it returns a -1
*/
int ReturnNameIndex(char **array, char *Str, int n)
{
	int i,idx,found=0;
	for(i=0;i<n;i++)  {
		if(strcmp(Str,array[i])==0) {
			if(found)
				return(-2);
			else {
				found =1;
				idx = i;
			}
		}
	}
	if(found) {
		return(idx);
	}
	
	return(-1);
}	