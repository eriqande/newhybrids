


/*  these are options for the type of simulation we are doing. */
#define NEW_HYBRIDS 0
#define PRITCHARD  	1
#define JABES		2
#define JABES_NO_PURES 3

/*  options for how to initialize the variables.  These are assigned to
the variable D->StartingDispersionTheta  */
#define PRIORS_ONLY 0
#define USE_PRIORIZED_ALLELES 1
#define USE_FIXED_ZS_4_PI_INIT 1


/* some global variables */
GLOB int gAlphasConstrainedEqual;  /*  1 for yes 0 for no. */

/* in case you want very strong prior on Pi---basically to fix it at any certain value...*/
GLOB int gPiFixed;  /* 1 for yes 0 for no */
GLOB double *gPiFixedValues;
GLOB char gPWD[10000];  /* global variable to store the PWD */


#define JA_PURE 0
#define JA_ADMIXED 1


#define MAX_ALLELES 999
#define MAX_GN      80000
#define MAX_LOCUS_NAME_LENGTH 80
#define MAX_CAT_NAME_LENGTH	  80
#define MAX_INDIV_NAME_LENGTH 80
#define MAX_DATAFILE_NAME_LENGTH 100
#define MAX_NUM_SOURCES_PRITCH  20  /*  maximum number of components (sub-pops) for the Pritchard method */
#define MAX_ALLELE_NAME_LENGTH 50
#define MAX_PATH_NAME_LENGTH 10000

#define CYCLE_alk(P) for(a=0;a<2;a++)  {for(l=0;l<P->L;l++)  {for(k=0;k<P->Kl[l];k++)  {
#define CYCLE_lka(P) for(l=0;l<P->L;l++)  {for(k=0;k<P->Kl[l];k++)  {for(a=0;a<2;a++)  {
#define CYCLE_ilc(P) for(i=0;i<P->M;i++)  {for(l=0;l<P->L;l++)  {for(c=0;c<2;c++)  {
#define END3CYCLE }}}

#define CYCLE_al(P) for(a=0;a<2;a++)  {for(l=0;l<P->L;l++)  {
#define DO_W_CYCLE(W) for(W[0]=0;W[0]<2;W[0]++) { for(W[1]=0;W[1]<2;W[1]++) {
#define CYCLE_TO_1(w) for(w=0;w<2;w++)  {
#define END2CYCLE  }}

#define CYCLE_i(P) for(i=0;i<P->M;i++)  {
#define CYCLE_g(P) for(g=0;g<P->Gn->v;g++)  {
#define CYCLE_l(P)   for(l=0;l<P->L;l++)  {
#define CYCLE_k(P) for(k=0;k<P->Kl[l];k++)  {
#define CYCLE_a for(a=0;a<2;a++)  {
#define CYCLE_c for(c=0;c<2;c++)  {
/*  here is one for cycling over the components in the Pritchard Methodm */
#define CYCLE_s(P) for(s=0;s<P->NumSources->v;s++)  { 
#define CYCLE_jv for(jv=0;jv<2;jv++)  {  /*  this one is for cycling over JA_PURE and JA_ADMIXED */


 
#define END1CYCLE  } 

/* these are Macros for the index of the individual option Flag */
#define NUM_INDIV_FLAGS 9
#define Z_FIXED 0
#define FIXED_Z 1
#define Z_FIXED_IN_PURE 2
#define PURE_POOL 3  /*  if Z_FIXED_IN_PURE == 1, this should be set to tell
						which of the gene pools, i.e. 0 or 1 this corresponds to */
#define DONT_CONTRIBUTE_TO_PI 4
#define COMPUTE_POFZ_ANYWAY 5
#define LOCALE 6 
#define HAS_LOCALE 7
#define DONT_CONTRIBUTE_TO_THETA 8


/*  an enum for different types of loci */
enum locus_type
{
	UNSET,			/*  Not yet set  */
	CODOM, 			/*  codominant */
	AFLP,			/*  */
	HAS_NULL		/*  a locus that has a null allele */
};

/*  an enum for different types of priors */
enum prior_type
{
	UNIFORM,		/*  0 is for a uniform Dirichlet prior */
	JEFFREYS,		/*  1 is for a Jeffreys Dirichlet prior */
	USER_DIRCH,		/*  2 is for a Dirichlet Prior with user-specified parameters */
	FIXED_PRIOR			/*  3 means that the value of the proportions are fixed (prior with infinite weight) */
};


/*  a struct to hold all the data, most of which won't get changed at all. */
typedef struct
{
	int *Kl;  	/*  the number of alleles in each locus */
	int L;		/*  the number of loci */
	int M;		/*  the number of individuals in the sample */
	
	int TypeOfSim;  /*  the declared type of simulation to be done at the beginning here */
	int StartingDispersionTheta;  /* for different ways of initializing the variables */
	int StartingDispersionPi;
	
	enum locus_type *LocTypes;  	/*  array of enums telling us what kinds of loci we have */
	
	int AFLPsPresent;  /* takes a 0 if no AFLP loci in the data set, 1 otherwise */

	char DataFileName[MAX_DATAFILE_NAME_LENGTH];
	char **LocNames;  /*  array of strings telling us the locus names */
	char ***AlleleNames;  	/*  array of array of strings giving the names of the alleles at all */
							/*  the loci. */
							
	ival ****Yobs;  	/*  the actual, OBSERVED data themselves, subscripted by [i][l][c]---i is the individual */
				/*  subscript, l is the locus subscript, c is the gene copy subscript (0 or 1, */
				/*  depending on if it is the first or second gene copy in the individual). */
				/*  values of Y which are >=0 will indicate different types of alleles, and */
				/*  == -1 will be our flag for a locus with missing data.   */
				/*  These OBSERVED values will be exactly what the Latent Data Y values are for */
				/*  codominant loci, but for other types of locus we shall have a difference */
				/*  FOR AFLP's Yobs[i][l][0] is 0 for  "-" individual and Yobs[i][l][0] is 1 for a "+" */
				/*  individual */
				
	ival *Gn;  /*  the number of genotype frequency classes  */
	
	char **CategoryNames;
	
	double ***G;	/*  the expected genotype frequencies of the Gn categories */
					/*  these two are included amongst the unkown in case I want */
					/*  to try sampling over them at some point.   */
					/*  G[g][0][0] is for 00, G[g][0][1] is for 01, */
					/*  G[g][1][0] is for 10, and G[g][1][1] is for 11 given the g-th category */
					/*  Note that the "/2" in equation 2 is already taken care of because we have */
					/*  G[g][0][1] = G[g][1][0].  ;   */
					
	int **IndFlags;
	
	char **IndNames;
	
	double ***PriorizedAlleleCounts;
	
	
	int PiTraceReport;
	int ZTraceReport;
} hyb_data;



/* I believe these are pretty well obsolete, being replaced by the #defined Flags */
typedef struct 
{
	int Z_Fixed;  /*  takes a 1 if the GtypFreqCategory is known, 0 otherwise.  Default is 0 */
	int Fixed_Z;	/*  takes the value that Z for the individual should be set to. */
	int FixedInPure;  /* takes a 1 if Fixed_Z is a Pure category, 0 otherwise */
	int ComputePofZAnyway;		/*  1 if the PofZ should be computed for the individual.  This should only */
							/*  come into play if SourceKnown == 1 or GtypFreqCatKnown == 1 but we */
							/*  still want to ComputePofZ to see how well we would have been able to */
							/*  place the individual given our data and assumptions. */
	int ComputeLatentYsForKnownSource;   /*   if the individual is of known GtypFreqCat but has */
										 /*  AFLP's and so needs the latent  Y's to get summed  */
										 /*  over then this should be set to 1 and that calculation */
										 /*  will be carried out for the Pure GtypFreqCat to which the  */
										 /*  individual belongs, but only at the AFLP loci. */
	int ContributeToPi;	/*  Under the Gtyp Freq Categories model, if this is 1, then the individual's Z */
						/*  gets counted in computing the full conditional for Pi.  If 0, then the  */
						/*  individual's Z doesn't count toward computing the full cond for Pi. Default is 1 */
	int AllelesContribute;  /*  1 if the alleles in the individual should be divvied up according */
							/*  to their W's when deriving the full conditional for the allele freqs. */
							/*  O if the individuals alleles should NOT contribute to the estimation  */
							/*  of frequencies.  I anticipate this may not really get used. */

}  hyb_ind_flags;  


/*  here is a simple struct for the latent variables for each individual that are */
/*  unique to Jonathan Pritchard's method. */
typedef struct 
{
	dval **Q;  /*  do this as a vector, so it extends to the multiple possible sources case easily. */
	ival **NumW;  /*  the number of genes, overall, in an individual, allocated to each component */
					/*  subscripted by [s] */
					
	dval **GtypProb;  /*  GtypProb[0] is given it is pure */
					  /*  GtypProb[1] is given it is admixed, and conditional on current */
					  /*  values of the alphas */
					  
	ival *V;  /*  the indicator of whether someting is pure (JA_PURE, 0) or admixed (JA_ADMIXED, 1) */
	dval **PofV;  /*  for storing the probabilities that an individual is pure or not. */
	
	
} Pritch_Indiv_Latent_struct;

typedef struct
{
	ival *NumSources;  /*  the number of components in the mixture */
	dval **Alpha;  /*  the alpha parameter---note that I will allow for non-symmetrical Dirichlet */
				  /*  priors on the Q variable */
	dval **MaxAlpha;  /*  max and min values that each component of alpha is allowed to take */
	dval **MinAlpha;
	
	dval **Sigma_SD;  /*   standard deviation of the normal proposal for elements of alpha.  This is a  */
					  /*  vector because we allow for different proposal widths for different components */

	dval **Xi;		/*  this is the mixing proportion between pures and admixed.  NOTE:  for Pi, we use the */
					/*  Pi that is already in the standard NewHybrids struct. */

	ival **Vcounts;  /*  the number of individuals currently in either the pure or admixed categories. */

} Pritch_Overall_Latent_struct;

/*  a struct to hold the latent variables that all apply to a single individual */
/*  in the sample */
typedef struct
{
	int idx;		/*  just stores the index of the individual.  This is useful for */
					/*  finding the Yobs for the individual in functions that just get */
					/*  a pointer to hyb_indiv passed to it.   */
					
	int *Flag;  /*  points to the flags for the individual that tell about its status. */
	
	ival ***Y;		/*  The unobserved data.  For codominant alleles, this will just be Yobs.  However, */
					/*  for other locus systems these values will change and be sampled over.  This is, */
					/*  however, the Y that all the probabilities of genotype will take place on. */
					/*  subscripted by [l][c] */
	
	ival ***W;	/*  an array similar in shape to that of Y in the data.  The elements indicate */
				/*  the species of origin of particular gene copies in the individuals */
				/*  subscripted by [l][c] where c = 0 or 1 */
	
	ival *Z;		/*  Z[i] points to an ival that indicates the genotype frequency */
					/*  class of the individual */
					
	dval **LogPofY;  /*  Log of the prob of the genotype given the current thetas and each of the different */
					/*  genotype frequency classes.  Subscripted by [g] */
					
	dval **PofZ;	/*  The probability that the individual belongs to each of the gtyp freq classes. */
					/*  subscripted by [g] */
					
	dval **UniPriPofZ;  /* this is a structure used to store the PofZ that would occur if the prior were
							uniform.  Basically this is a way of storing the likelihood that the individual
							belongs in one of the different hyb categories, and the likelihoods have been 
							scaled to sum to one.  The information in it is NOT used in making any further udates */
					
	dval ****PofW;  /*  For each locus [l], the probability that the W's are 00, 01, 10, or 11 */
					/*  	subscripted by [l][w1][w2] where w1 \in {0,1} and w2 \in {0,1}.  Currently, for */
					/*  the Pritchard case, this is not used */
					
	double *****PofWandYlat;  /*  for each locus [l] and origin of markers [w1] or [w2] and type */
							/*  of AFLP marker 0, or 1, at copies one and two of the locus, this */
							/*  will store the probability of that. */
							/*  Subscripted by [l][w1][w2][y1][y2] */
							
	/*  a pointer to the Pritchard method extra latent variables for individuals (basically just */
	/*  Q)  */
	Pritch_Indiv_Latent_struct *PritInd;
	
					
} hyb_indiv;


/*  a struct to hold all the latent data and unknown parameters */
/*  and also to hold some space that will get used repeatedly */
typedef struct 
{
	int TypeOfSim;  /*  for the type of simulation: */
					/*  0 NewHybrid  (My genetics paper) */
					/*  1 Pritchard method */
					/*  2 My method from JABES paper.  These are defined above as constants */
	
	ival *Gn; 		/*  the number of genotype frequency categories */
	
	double ***G;	/*  the expected genotype frequencies of the Gn categories */
					/*  these two are included amongst the unkown in case I want */
					/*  to try sampling over them at some point.   */
					/*  G[g][0][0] is for 00, G[g][0][1] is for 01, */
					/*  G[g][1][0] is for 10, and G[g][1][1] is for 11 given the g-th category */
					/*  Note that the "/2" in equation 2 is already taken care of because we have */
					/*  G[g][0][1] = G[g][1][0].   */
					
	dval **Pi;  	/*  the mixing proportions.  Subscripted by [g] */
	
	dval ****Theta;  /*  The allele frequencies in the different populations */
					 /*  subscripted [a][l][k] */
	
	double ***r;  	 /*  this is used along the way---it counts the number of alleles */
					 /*  of the different types currently allocated to each of the two spp. */
					 /*  and then the prior will be added to it.  Subscripted [a][l][k]. */
					 
	double *s;  		 /*  used to count up the number of individuals in each geno. freq. cat. */
					 /*  similar to r. */
					 
	
	hyb_indiv **Ind;  /*  pointer to an array of individuals	 */
	
	double TheCompleteDataLogLike;
	
	Pritch_Overall_Latent_struct *PritLat; 	
	
	dval **Locus_KB;   /*  an array of dvals subscripted by l.  Each carries the Kullback Leibler */
					   /*  divergence between the populations for each locus	 */
	
} hyb_latent;

typedef struct 
{
	dval ****Lambda;  	/*  the allele frequency priors */
	dval **Zeta;  		/*  the mixing proportion priors */
	
	dval **XiPrior;    /*  the priors on the mixing proportion of Pure vs Admixed in the JABES version */
	
	enum prior_type PiPriType;
	enum prior_type ThetaPriType;
} hyb_prior;

/*  A big monster struct to hold all info relative to a chain */
typedef struct
{
	hyb_data *Dat;
	hyb_latent *Lat;
	hyb_prior *Pri;
	
	int Seed1;
	int Seed2;

} hyb_chain;


/* to store options that can be passed on the command line */
typedef struct
{
	char DataFilePath[MAX_PATH_NAME_LENGTH];
	char GtypCatFilePath[MAX_PATH_NAME_LENGTH];
	double *GProbs[MAX_GN];  /* for the genotype probs of the categories GProbs[g][x] is the prob that 
						an indiv in category g has x copies of genes at a locus from Species 0 (x = 0,1,2).
						So, this is different than specified in the GtypsProb File, which I will do away with
						altogether.
						*/ 
	char *GProbNames[MAX_GN];
	int NumGProbs;
	char AlleFreqPriorPath[MAX_PATH_NAME_LENGTH];
	long Seed1;
	long Seed2;
	
	enum prior_type PiPriType;
	enum prior_type ThetaPriType;
	
	double *PiFixedValues;  /* to store a special user-set dirichlet prior */
	double *PiUserDirchValues;
	
	int NumBurnIn;
	int NumPostBurnIn;	
	
	int NoGui;
	
	int PiTraceReport;  /* The frequency with which the program should print out a trace of the Pi values.  0==None! and any number greater than 0 is the
							number of sweeps between reports.  */
    int ZTraceReport;
} cli_opts;



#define MAX_LOCALES  100
#define CLINE_POINTS 300
#define LOCATION_UNSET -999.887766 
/* some structure to hold information about clines at different loci */
typedef struct {
	dval **AlleFreq; /* observed frequencies at each locale at this locus */
	ival **Y;  /* the number of gene copies assigned to spp. 1 in each locale */
	ival **N;   /* the total number of gene copies assigned to either species in each locale */
	dval *alpha;  /* the cline center parameter */
	dval *beta;  /* the cline shape parameters */
	
	/* these are related to the metropolis scheme */
	dval *a_sd;
	dval *b_sd;
	dval *LogLike;  /* current value of the log-likelihood */
	
	
	/* these are to hold fitted values of the alle freqs (i.e. to plot the clines) */
	dval **ClineY;  /* for the current cline */
	
	
} locline_struct;

/* I will probably declare one of these as a single global variable for now */
typedef struct {
	int NumLocales;  /* the number of locales */
	int MaxLocaleIndex;  /* the value of the highest-indexed locale */
	int *HasLocale;  /* HasLocale[i]=1 if there is a locale with index i, 0 otherwise */
	double *x;  /* the position of the locales along the cline */
	int *NumInds; /* number of individuals in each locale */
	
	int NumClineXs;  /* number of points along the x-direction for plotting the clines */
	double *ClineX;  /* the x-values for plotting the clines.   */
	
	locline_struct **Locs;  /* for keeping all the locus-specific information */
	
} cline_struct;


/* here is our global variable */
GLOB cline_struct *gClines;



/*  The function prototypes */
extern void Count_r(hyb_chain *C);
extern void Count_s(hyb_chain *C);
extern void UpdateTheta(hyb_chain *C);
extern void UpdatePi(hyb_chain *C);

extern double StoreSingleLogPofY(hyb_indiv *Ind, hyb_chain *C);
extern double PofY_at_a_locus(hyb_indiv *Ind,
		hyb_chain *C, int l, int g);
extern void UpdateWandZandY(hyb_chain *C);
extern int IntegerFromProbsRV(int Num,dval **P);
extern void DrawNewW_and_LatentY_indiv(hyb_indiv *Ind, hyb_chain *C);
extern void SingleSweep(hyb_chain *C);
extern void InitializeChain(hyb_chain *C);
extern hyb_indiv *CreateIndiv(hyb_data *D, int i);
extern hyb_chain *CreateLatentChain(hyb_data *D, hyb_prior *P);
extern hyb_prior * CreatePriors(hyb_data *D, 
		enum prior_type PiPriorType, enum prior_type ThetaPriorType);

extern void SetPiPriors(hyb_data *D, hyb_prior *P, enum prior_type PiPriorType);
extern void SetThetaPriors(hyb_data *D, hyb_prior *P, enum prior_type LambdaPriorType);

extern void IncrementValues(hyb_chain *C, int DoThetas);
extern void OutputHistograms(hyb_chain *C);
extern void ResetAllAveragesEtc(hyb_chain *C);
extern int YnotMissing(hyb_chain *C, int i, int l);

extern double KullbLeibAtLocus(hyb_chain *C, int l);
extern void KullbLeib(hyb_chain *C);


/*  data stuff */
extern cli_opts *Get_NewHybs_CLI_Opts(int argc, char *argv[]);
extern int ReadIndivOptions(FILE *s, int *F, int Idx, char **Name);
extern void ProcessIndivOptions(hyb_data *D);
extern void PriorizeAllelesFromFixedZ(hyb_data *D);
extern void AddPriorizeAllelesFromFile(hyb_data *D, char *File);
extern int ReturnNameIndex(char **array, char *Str, int n);
extern hyb_data * GetData(char *FileName);
extern const char *LocType2str(enum locus_type loc);
extern const char *PriorType2str(enum prior_type loc);
extern void ReadDataLumpedFormat(int temp_int, int HasLocusNames, int ***TempY, 
				int TempM, int TempL, enum locus_type *TempLocTypes, int Divisor, FILE *in, int **IndivFlags, char **IndNames);
extern void ReadDataNonLumpedFormat(int temp_int, int HasLocusNames, int ***TempY, 
				int TempM, int TempL, enum locus_type *TempLocTypes, FILE *in, int **IndivFlags, char **IndNames);
extern void GetGtypFreqCats(hyb_data *D, char *FileName);
extern void CopyGtypFreqCatsFromCL(hyb_data *D, cli_opts *CL);
void SetUpDefaultGtypCategs(cli_opts *CL);


/*  Data output stuff */
extern void fprintCL_Probs(FILE *out, cli_opts *CL);
extern void fprint_PofZ(hyb_chain *C);
extern void fprint_AlleleAverages(hyb_chain *C);
extern void fprint_PiAverages(hyb_chain *C);
extern void fprint_UniPriPofZ(hyb_chain *C);
extern void printPi_Trace(int Rep, hyb_chain *C, int PrintHeader);
extern void printZ_Trace(int Rep, hyb_chain *C, int PrintHeader);

/* controls whether or not the GLUT windows pop up */
int RunWithoutGraphics(hyb_chain *C, int DoAsBurnIn, int DoAsReps);


/******************* 
CLINE RELATED STUFF
*************************/
void ProcessClineOptions(hyb_data *D, cline_struct *C);
locline_struct *AllocLocline(int MaxLocaleIndex,  double NumClineXs,  double asd,  double bsd);
void UpdateClineVariables(hyb_chain *C, cline_struct *L);
void RandoStartClineVars(hyb_chain *C, cline_struct *L);
void IncrementClineVars(hyb_chain *C, cline_struct *L);
void ResetClineVarAverages(hyb_chain *C, cline_struct *L);






/***************************************************************************
*
*		FUNCTIONS FOR PRITCHARD ET AL'S METHOD
*
***************************************************************************/
void PritchUpdateW(hyb_chain *C);
void PritchUpdateQ(hyb_chain *C);
Pritch_Indiv_Latent_struct *AllocPritchInd(int NumSources);
Pritch_Overall_Latent_struct *AllocPritchOverall(int NumSources);
void AllocatePritchardParts(hyb_chain *C);
void InitializePritchVars(hyb_chain *C);
void PritchUpdateAlphaRComp(hyb_chain *C, int r);
void PritchUpdateRandomCompOfAlpha(hyb_chain *C);
void IncrementPritchVars(hyb_chain *C, int DoThetas);
void PritchResetAllAveragesEtc(hyb_chain *C);
void PritchSingleSweep(hyb_chain *C);
void UpdateWbyBaum(hyb_chain *C);
void JABESComputePofZ(hyb_chain *C);
void JABESComputePofV(hyb_chain *C);
void JABESUpdateV(hyb_chain *C);
void JABEScount_Vcounts(hyb_chain *C);
void JABES_UpdateXi(hyb_chain *C);
void JABESUpdateAlphaUsingQ(hyb_chain *C, int r);
void JABESUpdateRandomCompOfAlpha(hyb_chain *C);
void JABES_UpdatePi(hyb_chain *C);
void JABES_UpdateZs(hyb_chain *C);
void JABES_UpdateWForPure(hyb_chain *C);
void JABES_StoreLogPofY_Pure(hyb_chain *C, int i, int a);
void UseLastFromNForPriorForZero(hyb_chain *C, int N);




















