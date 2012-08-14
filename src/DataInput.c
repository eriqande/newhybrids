#include <string.h>


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "ECA_MemAlloc.h"
#include "MCTypesEtc.h"
#include "ECA_utilities.h"
#include "ECA_ReadData.h"
#include "ECA_Opt3.h"
#include "MathStatRand.h"

#include "NewHybrids.h"



#define NH_LUMPED 0
#define NH_LINED 1
#define NH_NON_LUMPED 2


/*  prototypes for functions in this file: */
int NotInArray(int M, int *Array, int Num);
int ReturnSubscript(int *Array, int Length, int Compare);
void CheckLocType(char *inStr, enum locus_type *LocType, int i, int l);



/* this function is used to extract options off the command line.  It stores those options
in a cli_opts struct, from whence they may be copied over to the program data structures */
cli_opts *Get_NewHybs_CLI_Opts(int argc, char *argv[])
{
	cli_opts *temp = (cli_opts *)ECA_MALLOC(sizeof(cli_opts));
	int datpathF=0,
	gtyp_cat_fileF=0,
	gprobsF=0,
	afreq_priorF=0,
	theta_priorF=0,
	pi_priorF=0,
	piprivec_f=0,
	seedF=0,
	numburninF=0,
	numpostburnF=0,
	noguiF=0,
	NoMoreCategories=0,
	print_tracesF=0;
	
	
	DECLARE_ECA_OPT_VARS;
	
	
	
	SET_OPT_WIDTH(24);
	SET_ARG_WIDTH(16);
	SET_PROGRAM_NAME("newhybrids");
	SET_PROGRAM_SHORT_DESCRIPTION("Computes posterior probability of hybrid categories");
	SET_VERSION("VERSION: 2.0+ Developmental.  July/August 2007");
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson (eric.anderson@noaa.gov)");
	SET_VERSION_HISTORY("");
	SET_PROGRAM_LONG_DESCRIPTION(\tNewHybrids is a program for computing the posterior distribution that individuals in a 
		sample fall into different hybrid categories. The details of the algorithm are published in
		Anderson and Thompson \0502002\051 Genetics 160:1217-1229.\n\n\tThere are two flavors of the program.
		The plain vanilla version has no graphical user interface\054 and runs with a simple command line
		interface.  The fancier version\054 though still launched off the command line\054 
		has a graphical user interface based on OpenGL and GLUT libraries. This 
		graphical interface  allows the real-time observation of most of the variables involved in 
		the Markov chain Monte Carlo \050MCMC\051 simulation. This is valuable for assessing the relia- 
		bility of the results from the MCMC run.  In order to take full advantage of these features\054 
		you will want to read "A Users Guide to the GLUT for Markov Chain Monte Carlo Graphical 
		Interface" which is included in the standard distribution of the NewHybrids software as the file 
		GFMCMC_UserGuide.pdf.  Disabling the GLUT interface is done with the --no-gui option.);

		
	
	/* set some default flags for some of the variables */
	temp->GtypCatFilePath[0] = '\0';  /* set it to a zero-length string */
	temp->NumGProbs = 0;
	temp->PiFixedValues=NULL;
	temp->PiUserDirchValues=NULL;
	temp->NumBurnIn = 10000;
	temp->NumPostBurnIn = 50000;
	sprintf(temp->DataFilePath,"UNSET");
	sprintf(temp->AlleFreqPriorPath,"UNSET");
	temp->PiPriType=JEFFREYS;
	temp->ThetaPriType=JEFFREYS;
	temp->NoGui = 0;
	temp->PiTraceReport = 0;


	BEGIN_OPT_LOOP
	
	
		OPEN_SUBSET(
			Command Line Switches for Standard Analysis,
			Standard analysis settings,
			These standard settings control the data used and how NewHybrids analyzes it. 
		);
		if (  REQUIRED_OPTION(
			Data File,
			datpathF,
			d,
			data-file,
			F,
			pathname to the data file,
			The name of the data file that you want NewHybrids to analyze.  You can specify it as a full pathname.
			The format of the data file is described in the NewHybrids documentation.  
			) )	{
			if(ARGS_EQ(1)) {
				GET_STR(temp->DataFilePath);
			}
		}
		if(OPTION(
			Gtyp Categ File,
			gtyp_cat_fileF,
			c,
			gtyp-cat-file,
			F,
			path to file holding the genotype category probabilities ,
			Use this to specify a file holding genotype category probabilities as described in the manual. 
			The standard set to use is all possible products of two generations of mating---that is\054 Pure_A\054
			Pure_A\054 F1\054 F2\054 BX_A\054 BX_B.  If this option is omitted\054 and the -g/--gtyp-freq-probs is not
			given\054 then that standard analysis will
			be performed.  If both this option and the -g/--gtyp-freq-probs option are given\054 the program will quit
			with an error.)) {
			
				if(gprobsF>0) {
					fprintf(stderr,"Error! You issued both the -c/--gtyp-cat-file option AND a -g/--gtyp-freq-probs option.  They cannot both be issued.  Exiting\n");
					exit(1);
				}
				if(ARGS_EQ(1)){  char GtypCatFile[10000];
					GET_STR(temp->GtypCatFilePath);
					
				}
			
		}
		if( MULT_USE_OPTION(
				CL-Specified Genotype Probs,
				gprobsF, 
				g,
				gtyp-freq-probs, 
				S R0 R1 R2,
				specify genotype frequency category S with probabilities R0 R1 and R2,
				Use this option to define ON THE COMMAND LINE the hybrid categories individuals may belong to. 
				This may or MAY NOT be more convenient than specifying them in a file called with the -c/--gtyp-cat-file
				option.
				Assume there may be hybrids between a species A and a speceis B. 
				S is the name given to the category.  R0 is the expected proportion of loci drawn from
				a member of the category having 0 genes from species B.  R1 is the proportion of loci having
				exactly one copy from species B.  And R2 is the proportion of loci having 2 copies from species B.
				For example: "-g Pure_A 1.0 0.0 0.0" would define the category for pure species A.  Alternatively
				"-g F1 0.0 1.0 0.0" would define the category for F1.  "-g BX_A .5 .5 0" defines the category of a 
				backcross to the A species.  You should follow the convention that the first two categories you define are
				the pure A and pure B categories.  The progam checks for categories that are duplicates---either having the same name
				or near the same expected proportions.  Note that the format for the arguments to this option is DIFFERENT than the rows
				in the file called by --gtyp-cat-file. If both this option and the -c/--gtyp-cat-file option are given\054 the program will quit
			with an error. ,
				MAX_GN
				) ) { 
				if(gtyp_cat_fileF>0) {
					fprintf(stderr,"Error! You issued both the -c/--gtyp-cat-file option AND a -g/--gtyp-freq-probs option.  They cannot both be issued.  Exiting\n");
					exit(1);
				}
				if(ARGS_EQ(4)) {
					int i;
					double normo;
					
					if(NoMoreCategories) {
						fprintf(stderr,"Error! You have already listed prior choices for hybrid categories using the \"Dirichlet\" or the \"fixed\" keywords with option --pi-prior.  You cannot add more categories using the -g/--gtyp-freq-probs option after this.  Exiting...\n");
						exit(1);
					}
					
					/* first get the name of the category */
					temp->GProbNames[temp->NumGProbs] = (char *)ECA_CALLOC(MAX_CAT_NAME_LENGTH,sizeof(char));
					GET_STR(temp->GProbNames[temp->NumGProbs]);
					
					/* check for duplicate gtyp freq category names */  
					for(i=0;i<temp->NumGProbs;i++)  {
						if(strcmp(temp->GProbNames[i],temp->GProbNames[temp->NumGProbs])==0) {
							fprintf(stderr,"Error!  The name %s for a genotype frequency category has been used more than once.\n",temp->GProbNames[temp->NumGProbs]);
							OPT_ERROR;
						}
					}
					
					/* now, collect the genotype freq category probs */
					temp->GProbs[temp->NumGProbs] = (double *)ECA_CALLOC(3,sizeof(double));
					for(normo=0.0,i=0;i<3;i++)  {
						temp->GProbs[temp->NumGProbs][i] = GET_DUB;
						normo += temp->GProbs[temp->NumGProbs][i];
					}
					/* now normalize them */
					for(i=0;i<3;i++)  {
						temp->GProbs[temp->NumGProbs][i] /= normo;
					}
					
					/* now check to see if it is like any others */
					for(i=0;i<temp->NumGProbs;i++)  { int j; double diff;
						diff=0;
						for(j=0;j<3;j++)  {
							diff += ECA_ABS(temp->GProbs[temp->NumGProbs][j] - temp->GProbs[i][j]);
						}
						if(diff<.01)  {
							fprintf(stderr,"Error! Category %s genotype freqs are almost identical to those of category %s.  They won't be estimable.  This is considered an error.\n",temp->GProbNames[i],temp->GProbNames[temp->NumGProbs] );
							OPT_ERROR;
						}
					}
					
					/* increment the number of categories counter */
					temp->NumGProbs++;
					
				}
		}
		if(OPTION(
			Allele Frequency Priors File,
			afreq_priorF,
			,
			alle-prior-file,
			F,
			path to the file holding information about the priors for allele frequencies,
			F is the name of the file holding the priors for allele frequencies.  The format of that file is
			described in the NewHybrids documentation.  To form the file it is best to run NewHybrids on the data set
			using the default prior and then modify the resulting file aa-LociAndAlleles.txt renaming it to whatever
			you wish to call it.  You supply its filename as the string F to this option.  This option is not required.
			If it is omitted then inference of the allele frequencies depends only on the the --theta-prior option.))  
		{
			if(ARGS_EQ(1)) {
				GET_STR(temp->AlleFreqPriorPath);
			}
		}
		if(OPTION(
			Allele Frequency Prior Type,
			theta_priorF,
			,
			theta-prior,
			C{uniform;Jeffreys},
			Set the type of prior used for the prior on theta,
			This option sets the type of prior used for theta---the allele frequencies.
			S is a string that denotes the choice.  Currently there are only two options: "uniform"
			and "Jeffreys".  Therefore S should be either "uniform" or "Jeffreys". 
			Both choices are Dirichlet prior distributions.  The uniform prior is a Dirichlet distribution with
			all elements of the parameter vector equal to 1.  The Jeffreys prior is the same except each element of
			the Dirichlet parameter vector is 1 divided by the number of alleles at the locus.  The default is Jeffreys. 
			))
		{
			if(ARGS_EQ(1)) { char str[500];
				GET_STR(str);
				if(strcmp(str,"uniform")==0) {
					temp->ThetaPriType = UNIFORM;
				}
				else if(strcmp(str,"Jeffreys")==0) {
					temp->ThetaPriType = JEFFREYS;
				}
				else {
					fprintf(stderr,"Error! String argument \"%s\" to option --theta-prior not recognized.  S must be \"uniform\" or \"Jeffreys\". Exiting...\n",str);
					exit(1);
				}
			}
		}

		if(OPTION(
			Mixing Proportion Prior,
			pi_priorF,
			,
			pi-prior,
			C{uniform;Jeffreys;Dirichlet;fixed},
			Set the type of prior used for the prior on pi,
			Several different types of prior are available for pi---the mixing proportions
			of the different categories.  They are "uniform"; "Jeffreys"; "Dirichlet"; and "fixed".
			The uniform and the Jeffreys priors are special cases of the Dirichlet prior. 
			The uniform prior is a Dirichlet distribution with
			all elements of the parameter vector equal to 1.  The Jeffreys prior is the same except each element of
			the Dirichlet parameter vector is 1 divided by the number of hybrid categories.  If S=Dirichlet then 
			you have to supply the values of the Dirichlet parameter vector that you wish to use.  These values are R1
			up to Rk each one corresponding to one of the k hybrid categories that you must define using 
			k instances of the -g/--gtyp-freq-probs option BEFORE issuing the pi-prior option.  Using the S=fixed option also
			requires that the hybrid categories be specified earlier on the command line using the
			-g/--gtyp-freq-probs option.  S=fixed must also be followed by k arguments R1 up to Rk.  In this case the
			proportion of individuals in the population belonging to the j-th hybrid category will be fixed at Rj/[R1+...+Rk].

			))
		{	char str[500];
			
			GET_STR(str);
			if(strcmp(str,"uniform")==0) {
				temp->PiPriType = UNIFORM;
				if(COUNT_ARGS>0) {
					fprintf(stderr,"Error! There seem to be arguments remaining to in the --pi-prior option after \"uniform\".  There should not be.  Exiting...\n");
					exit(1);
				}	
			}
			else if(strcmp(str,"Jeffreys")==0) {
				temp->PiPriType = JEFFREYS;
				if(COUNT_ARGS>0) {
					fprintf(stderr,"Error! There seem to be arguments remaining to in the --pi-prior option after \"Jeffreys\".  There should not be.  Exiting...\n");
					exit(1);
				}	
			}
			else if(strcmp(str,"Dirichlet")==0) {
				temp->PiPriType = USER_DIRCH;
				if(ALREADY_HAS(gprobsF, "-g/gtyp-freq-probs")) {
					if(COUNT_ARGS == temp->NumGProbs) { int i;
						temp->PiUserDirchValues = (double *)ECA_CALLOC(temp->NumGProbs,sizeof(double));
						for(i=0;i<temp->NumGProbs;i++)  {
							temp->PiUserDirchValues[i] = GET_DUB;
						}
						NoMoreCategories = 1;
					}	
					else {
						fprintf(stderr,"Error! Incorrect number of arguments R1 ... Rk after keyword \"Dirichlet\" in option --pi-prior.  There should be %d.  Exiting...\n",temp->NumGProbs);
						exit(1);
					}	
				}
			}
			else if(strcmp(str,"fixed")==0) {
				temp->PiPriType = FIXED_PRIOR;
				if(ALREADY_HAS(gprobsF, "-g/gtyp-freq-probs")) {
					if(COUNT_ARGS == temp->NumGProbs) { int i; double normo;
						temp->PiFixedValues = (double *)ECA_CALLOC(temp->NumGProbs,sizeof(double));
						for(normo=0.0,i=0;i<temp->NumGProbs;i++)  {
							temp->PiFixedValues[i] = GET_DUB;
							normo += temp->PiFixedValues[i];
						}
						for(i=0;i<temp->NumGProbs;i++)  {
							temp->PiFixedValues[i] /= normo;
						}

						NoMoreCategories = 1;
					}	
					else {
						fprintf(stderr,"Error! Incorrect number of arguments R1 ... Rk after keyword \"fixed\" in option --pi-prior.  There should be %d  Exiting...\n",temp->NumGProbs);
						exit(1);
					}	
				}
			}
			else {
				fprintf(stderr,"Error! String argument \"%s\" to option --pi-prior not recognized.  S must be \"uniform\", \"Jeffreys\", \"Dirichlet\", or \"fixed\". Exiting...\n",str);
				exit(1);
			}
		}
		if(OPTION(
			Pi Prior Param Vector,
			piprivec_f,
			,
			pi-prior-vec,
			R1...Rk,
			Parameters for pi prior if choosing Dirichlet or fixed with option pi-prior ,
			If you chose Dirichlet or fixed with the --pi-prior option\054  then you must provide a vector of parameters
			for that prior.   These values are R1
			up to Rk\054 each one corresponding to one of the k hybrid categories that you must define using 
			k instances of the -g/--gtyp-freq-probs option or by using the --gtyp-freq-file option\054 either of which 
			must be issued BEFORE issuing the --pi-prior-vec option.  If you chose the Dirichlet prior\054 then 
			R1...Rk are the parameters of that Dirichlet prior distribution on pi.  If you chose the fixed prior\054 
			then the proportion of individuals in the population belonging to the j-th hybrid category will be fixed at Rj/[R1+...+Rk]. )) {
				if(ALREADY_HAS(pi_priorF,"--pi-prior") && ALREADY_HAS(gprobsF, "-g/gtyp-freq-probs")) {
					if(temp->PiPriType == USER_DIRCH) {
						if(COUNT_ARGS == temp->NumGProbs) { int i;
							temp->PiUserDirchValues = (double *)ECA_CALLOC(temp->NumGProbs,sizeof(double));
							for(i=0;i<temp->NumGProbs;i++)  {
								temp->PiUserDirchValues[i] = GET_DUB;
							}
							NoMoreCategories = 1;
						}	
						else {
							fprintf(stderr,"Error! Incorrect number of arguments R1 ... Rk in option --pi-prior.  There should be %d.  Exiting...\n",temp->NumGProbs);
							exit(1);
						}	
					}
					else if(temp->PiPriType == FIXED_PRIOR) {
						if(COUNT_ARGS == temp->NumGProbs) { int i; double normo;
							temp->PiFixedValues = (double *)ECA_CALLOC(temp->NumGProbs,sizeof(double));
							for(normo=0.0,i=0;i<temp->NumGProbs;i++)  {
								temp->PiFixedValues[i] = GET_DUB;
								normo += temp->PiFixedValues[i];
							}
							for(i=0;i<temp->NumGProbs;i++)  {
								temp->PiFixedValues[i] /= normo;
							}
							NoMoreCategories = 1;
						}	
						else {
							fprintf(stderr,"Error! Incorrect number of arguments R1 ... Rk in option --pi-prior-vec.  There should be %d  Exiting...\n",temp->NumGProbs);
							exit(1);
						}	

					}
					else {
						fprintf(stderr,"Error! You issued the --pi-prior-vec option but did not choose either Dirichlet or fixed pi priors!  Exiting...\n");
						exit(1);
					}
				}
			}
			
		if(OPTION(
			MCMC Burn In Length,
			numburninF,
			,
			burn-in,
			J,
			run for J sweeps of burn-in,
			The program will run for J sweeps and then reset the averages that it is collecting.  This serves to 
			allow the chain to "dememorize" its initial state.  If you are running NewHybrids with the graphical 
			output using the GLUT interface then this option is ignored and the averages must be reset manually.
			This option is optional.  The default value is 10000.
			))
		{
			if(ARGS_EQ(1)) {
				temp->NumBurnIn = GET_INT;
			}
		}
		if(OPTION(
			MCMC Run Length,
			numpostburnF,
			,
			num-sweeps,
			J,
			run for J sweeps AFTER burn-in,
			After burn-in the program will run for J sweeps of sample collection before terminating.  During this time
			the computed Monte Carlo averages in the output files are updated every 300 sweeps.  If the program terminates
			early the output to that point should still be in the files.  If you are running NewHybrids with the graphical 
			output using the GLUT interface then this option is ignored and the program will continue running until terminated 
			by the user.
			This option is optional.  The default value is 50000
			))
		{
			if(ARGS_EQ(1)) {
				temp->NumPostBurnIn = GET_INT;
			}
		}
		if(OPTION(
			Random Number Seeds,
			seedF,
			s,
			seeds,
			J1 J2,
			seeds for the random number generator,
			J1 and J2 must be positive integers.  If this option is not given the program
			will look for seeds in the file "newhyb_seeds" in the current working directory.
			If that file is not found then seeds will be generated by using the current time. 
			))
		{
			if(ARGS_EQ(2)) {
				temp->Seed1 = (long)GET_INT;
				temp->Seed2 = (long)GET_INT;
				if(temp->Seed1<=0 || temp->Seed2<=0) {
					fprintf(stderr,"Error! The two arguments to the -s/--seeds option must be greater than zero.  They can't be %ld and %ld.  Exiting...\n",temp->Seed1,temp->Seed2);
				}
				else {
					
				}
			}
			
		}
		if(OPTION(
			Disable GUI,
			noguiF,
			,
			no-gui,
			,
			disable the GLUT/OpenGL MCMC visualizer,
			Issuing this option disables the GLUT/OpenGL MCMC visualizer.  Instead output will be 
			what you get from the Without Grapics version.  Note that if you do not enter this option
			then you will enter the GLUT GUI and the program will run until you manually quit it.  You also
			have to manually reset the averages---there is no automatic burn in period.  See the GFMCMC guide
			for more information. )) {
				if(ARGS_EQ(0)) {
					temp->NoGui = 1;
				}
		}
		CLOSE_SUBSET;
	
		OPEN_SUBSET(Controlling Output,Controllin Output, Control Output of Traces)
		if(MULT_USE_OPTION(
						   Print Traces,
						   print_tracesF,
						   ,
						   print-traces,
						   S J ... ,
						   Tell program to print trace of variable type S every J sweeps,
						   Tell program to print trace of variable type S every J sweeps.  Currently S can only be the string Pi
						   and J tells how often to print the trace of the Pi variable.  J=1 means every sweep. J=5 means every 
						   five sweeps.  J must be 1 or greater. Later uses of this option overwrite previously set values., 10000)) {
			if(ARGS_GEQ(2)) { char temps[1000];
				GET_STR(temps);
				if( strcmp(temps,"Pi")==0 ) {
					if(COUNT_ARGS==1) {
						temp->PiTraceReport = GET_INT;
						if(temp->PiTraceReport < 1)  {
							fprintf(stderr,"Error! J option after Pi in --print-traces must be >=1, not %d. Exiting\n",temp->PiTraceReport);
							exit(1);
						}
					}
					else {
						fprintf(stderr,"Error! Expecting a single integer argument after the \"Pi\" argument to option --print-traces.\n",temps);
						OPT_ERROR;
					}
				}
				else {
					fprintf(stderr,"Error! Unrecongized S argument \"%s\" to option --print-traces.  Expecting the string \"Pi\"\n",temps);
					OPT_ERROR;
				}
			}
		}
		CLOSE_SUBSET
				  		
	END_OPT_LOOP
	
	/* down here, if we didn't get a gtyp_cat_file and we didn't get any -g options, then we set up the default, two-generation
	GProbs*/
	if(gtyp_cat_fileF + gprobsF == 0) {
		SetUpDefaultGtypCategs(temp);
	}
		
	return(temp);

}


/* this just populates the command line struct with the 6 "default" two generation genotype categories 
and their associated probabilities */
void SetUpDefaultGtypCategs(cli_opts *temp) 
{
	int i;
	
	temp->NumGProbs = 6;
	
	for(i=0;i<6;i++)  {
		temp->GProbNames[i] = (char *)ECA_CALLOC(MAX_CAT_NAME_LENGTH,sizeof(char));
		temp->GProbs[i] = (double *)ECA_CALLOC(3,sizeof(double));
	}
	sprintf(temp->GProbNames[0],"Pure_0");
	temp->GProbs[0][0]=1.0;   temp->GProbs[0][1]=0.0;  temp->GProbs[0][2]=0.0;
	
	sprintf(temp->GProbNames[1],"Pure_1");
	temp->GProbs[1][0]=0.0;   temp->GProbs[1][1]=0.0;  temp->GProbs[1][2]=1.0;
	
	sprintf(temp->GProbNames[2],"F1");
	temp->GProbs[2][0]=0.0;   temp->GProbs[2][1]=1.0;  temp->GProbs[2][2]=0.0;
	
	sprintf(temp->GProbNames[3],"F2");
	temp->GProbs[3][0]=0.25;   temp->GProbs[3][1]=0.5;  temp->GProbs[3][2]=0.25;
	
	sprintf(temp->GProbNames[4],"0_BX");
	temp->GProbs[4][0]=0.5;   temp->GProbs[4][1]=0.5;  temp->GProbs[4][2]=0.0;
	
	sprintf(temp->GProbNames[5],"1_BX");
	temp->GProbs[5][0]=0.0;   temp->GProbs[5][1]=0.5;  temp->GProbs[5][2]=0.5;
	
}

/* function to be called at or near the beginning of a line of data for
	an individual that will determine its options with respect to whether
	it has fixed z or not, etc.  Takes the input stream s, and puts the options
	in the indiv struct 
*/
int ReadIndivOptions(FILE *s, int *F, int Idx, char **Name)
{
	char p;
	int num=0,i;
	char pa[10];
	int ReadAPlace=0;
	
	/* read the next non-whitespace character.  If it is an integer or + or -, then it was a locus designation.  Put it back
	   on the stream to be read, and get out of this function, returning 0 (meaning no options were processed */
	for(p=getc(s);isspace(p);p=getc(s)) ; /* p is now the first non-whitespace character */
	if(isdigit(p) || p == '+' || p == '-') {  
		ungetc(p,s);  
		return(0);
	}
	/* then process the options, until it eats an int or a + or - when it could have gotten an option s, z or y */
	do {
		switch(p)  {
			case('z'):  /* Z_fixed option.  The fixed z value follows the z.  It must be followed by whitespace */
				F[Z_FIXED] = 1;
				for(p=getc(s);isspace(p);p=getc(s)) ;
				if(!isdigit(p))  { /* if the next character is not an integer given an error message */
					fprintf(stderr,"\nFatal error!  Option z used, but not followed by an integer");
					fprintf(stderr,"\nspecifying the gtyp freq category, at Individual %d.  Exiting...\n\n", Idx+1);
					exit(1);
				}
				else {   /* read the digits following z into a string and convert to a number */
					pa[0] = p;
					for(i=0;isdigit(pa[i]);i++)  {
						pa[i+1] = getc(s);
					}
					ungetc(pa[i],s);
					pa[i] = '\0';
					F[FIXED_Z] = atoi(pa);
				}
				break;
			case('s'):
				F[DONT_CONTRIBUTE_TO_PI] = 1;
				break;
			case('t'):
				F[DONT_CONTRIBUTE_TO_THETA] = 1;
				break;
			case('y'):
				F[COMPUTE_POFZ_ANYWAY] = 1;
				break;
			case('n'):  /* this is to give the individual a name */
				for(p=getc(s);isspace(p);p=getc(s)) ; /* eat any whitespace that might be there */
				(*Name) = (char *)ECA_CALLOC(MAX_INDIV_NAME_LENGTH,sizeof(char));  /* allocate space to it (thus making it non-NULL) */
				(*Name)[0] = p;
				for(i=0;!isspace((*Name)[i]);i++) {  /* cycle until you find a space */
					(*Name)[i+1] = getc(s);
				}
				ungetc((*Name)[i],s);
				(*Name)[i] = '\0';
/*	printf("\n\nGOT THE NAME:  %s\n",(*Name)); */
				break;
			case('p'):  /* The Place option (locale).  Place index must follow this.  It must be followed by whitespace */
				for(p=getc(s);isspace(p);p=getc(s)) ; /* eat any whitespace that might be there */
				if(!isdigit(p))  { /* if the next character is not an integer give an error message */
					fprintf(stderr,"\nFatal error!  Option p used, but not followed by a non-negative integer");
					fprintf(stderr,"\nspecifying the locale, at Individual %d.  Exiting...\n\n", Idx+1);
					exit(1);
				}
				else {   /* read the digits following p into a string and convert to a number */
					pa[0] = p;
					for(i=0;isdigit(pa[i]);i++)  {
						pa[i+1] = getc(s);
					}
					ungetc(pa[i],s);
					pa[i] = '\0';
					F[LOCALE] = atoi(pa);
					F[HAS_LOCALE] = 1;
					
					/* now, some bookkeeping */
					gClines->NumInds[F[LOCALE]]++;
					ReadAPlace = 1;
					if(gClines->HasLocale[F[LOCALE]]==0) {
						gClines->HasLocale[F[LOCALE]] = 1;
						gClines->NumLocales++;
						if(F[LOCALE] > gClines->MaxLocaleIndex) {
							gClines->MaxLocaleIndex = F[LOCALE];
						}
					}
				}
				break;
			case('x'):  /* Locale x-position option.  Place location must follow this.  It must be followed by whitespace */
				if(ReadAPlace == 0) {
					fprintf(stderr,"Reached the x flag for locale location without first getting the p flag\n");
					fprintf(stderr,"\nspecifying the locale position, at Individual %d.  Exiting...\n\n", Idx+1);
				}
				for(p=getc(s);isspace(p);p=getc(s)) ; /* eat any whitespace that might be there */
				if(!(isdigit(p) || p=='-' || p=='.') )  { /* if the next character is not a digit, a - or a . then bail*/
					fprintf(stderr,"\nFatal error!  Option x used, but apparently not followed by a real number");
					fprintf(stderr,"\nspecifying the locale position, at Individual %d.  Exiting...\n\n", Idx+1);
					exit(1);
				}
				else {   /* read the digits following p into a string and convert to a number */
					pa[0] = p;
					for(i=0;(isdigit(pa[i]) || pa[i]=='.' );i++)  {
						pa[i+1] = getc(s);
					}
					ungetc(pa[i],s);
					pa[i] = '\0';
					gClines->x[F[LOCALE]] = atof(pa);
					
				}
				break;

			default:
				break;
		}
		num++;
		for(p=getc(s);isspace(p);p=getc(s)) ;
	} while(isdigit(p)==0 && p != '+' && p != '-');
	
	/* put that last character back on the stream to be read as a locus, and then return the number of options processed */
	ungetc(p,s);
	return(num); 
}


/*  a function that looks at the IndFlags that have been read and then sets a couple of other
	ones as appropriate based on the values that were read */  
void ProcessIndivOptions(hyb_data *D)
{
	int i,g, PURE_0 = -1, PURE_1 = -1;
	
	/* first thing we have to do is figure out which of the categories is Pure  */
	CYCLE_g(D)
		if(D->G[g][0][0] == 1.0)
			PURE_0 = g;
		if(D->G[g][1][1] == 1.0)
			PURE_1 = g;
	END1CYCLE
	
	/* then cycle over the individuals and set their Z_FIXED_IN_PURE and PURE_POOL flags */
	CYCLE_i(D)
		D->IndFlags[i][PURE_POOL] = -1;  /* PURE_POOL should be -1 (meaningless) unless it gets set to a new value */
		if(D->IndFlags[i][Z_FIXED]) {
			if(D->IndFlags[i][FIXED_Z]==PURE_0)  {
				D->IndFlags[i][Z_FIXED_IN_PURE] = 1;
				D->IndFlags[i][PURE_POOL] = 0;
			}
			if(D->IndFlags[i][FIXED_Z]==PURE_1)  {
				D->IndFlags[i][Z_FIXED_IN_PURE] = 1;
				D->IndFlags[i][PURE_POOL] = 1;
			}
		}
	END1CYCLE
}

int NotInArray(int M, int *Array, int Num) 
{
/*
	Checks to see if an integer is in an array of ints.
	Returns a 1 if M is not in the array.  0 otherwise.
*/
	int i;
	int NotInThere = 1;
	
	for(i=0;i<Num;i++)  {
		if(M == Array[i]) {
			NotInThere = 0;
			break;
		}
	}
	return(NotInThere);
}

/*  return the position of the first occurrence of Compare in Array of length Length */
/*  Returns a -1 if Compare is not in the array. */
int ReturnSubscript(int *Array, int Length, int Compare)
{
	int i;
	int temp=-1;
	
	for(i=0;i<Length;i++)  {
		if(Array[i] == Compare)	{
			temp = i;
			break;
		}
	}
	
	return(temp);
}

/*  interprets the locus input string and decides what locus type it is. */
/*  If the LocType is UNSET, it sets it to the appropriate type.  If the  */
/*  LocType is set to a different type, it barks a warning */
void CheckLocType(char *inStr, enum locus_type *LocType, int i, int l)
{	
	enum locus_type TempLT = UNSET;  	/*  start with the assumption that it is UNSET */
	
	if(strcmp(inStr,"+")==0 || strcmp(inStr,"-")==0) {  /* this means the current locus is definitely AFLP */
		TempLT = AFLP;
	}
	else if(strcmp(inStr,"0")!=0) {  /* if not AFLP and not 0, then it is CODOM */
		TempLT = CODOM;
	}
	
	if(TempLT != *LocType)  {  /*  action of some sort will be necessary */
		if(*LocType == UNSET) {  /*  if the locus itself is currently unset, then we can set it */
			*LocType = TempLT;
		}
		else if(TempLT != UNSET){  /*  if the locus was not UNSET but TempLT != *LocType, then we report an error */
			 printf("\n\nError!  Locus %d, in Indiv %d, was %s.  But we just read",
			 		l,i,LocType2str(*LocType));
			 printf("a \'%s\' suggesting the locus is of type %s\n\n",inStr,LocType2str(TempLT));
		}
	}
}


hyb_data * GetData(char *FileName)
{
	int i,l,c,temp_int,k;
	int ***TempY;  /*  array of alleles carried at all loci by all individuals */
	int TempL,TempM, TempPolyL;	   /*  temporary number of loci and individuals */
	int *TempKl;  /*  temporary number of alleles at each locus */
	FILE *in,*out;
	int **TempIdx;
	hyb_data *D;
	int locus_counter;
	char tempStr[500];  /*  a temporary variable to hold strings */
	char **TempLocNames;  /*  temporary, to hold locus names */
	enum locus_type *TempLocTypes;
	int HasLocusNames = 0;  /*  flag for whether it has locus names or not */
	int Digits, Divisor, matched;
	char Format[100];  /*  to record the format of the data */
	int TheFormat;
	int **TempFlags;
	char ***TempAlleNames;
	char **TempIndNames;
	
	
	/* some Clines memory allocation */
	gClines = (cline_struct *)ECA_MALLOC(sizeof(cline_struct));
	gClines->HasLocale = (int *)ECA_CALLOC(MAX_LOCALES,sizeof(int));
	gClines->NumInds = (int *)ECA_CALLOC(MAX_LOCALES,sizeof(int));
	gClines->x = (double *)ECA_CALLOC(MAX_LOCALES,sizeof(double));
	gClines->NumClineXs = CLINE_POINTS;
	gClines->ClineX = (double *)ECA_CALLOC(gClines->NumClineXs,sizeof(double));
	gClines->NumLocales = 0;
	gClines->MaxLocaleIndex=0;
	for(i=0;i<MAX_LOCALES;i++)  {
		gClines->x[i] = LOCATION_UNSET;
	}	
	
	
	/*  first get all the data into the Temp structures */
	in = erdOpenFileOrRetry(FileName, "r");
	
	
	if( (matched = fscanf(in,"NumIndivs %d  NumLoci %d Digits  %d  Format %s",
				&TempM,&TempL, &Digits, Format)) != 4) {
		printf("\n\nSorry, file %s does not begin in the format:",FileName);
		printf("\n\nNumIndivs XX\nNumLoci XX\nDigits XX\nFormat \"Format Type\" (i.e. Lumped of Lined)\n");
		printf("\n\nGo back and fix that and try again");
		printf("\n\nExiting to system\n\n");
		exit(1);
	}
	
	/*  set the value of the divisor */
	Divisor = (int)(pow(10.0f,Digits) + .001f);
	
	/*  set TheFormat */
	if( strcmp(Format,"Lumped")==0  || strcmp(Format,"LUMPED")==0  || strcmp(Format,"lumped")==0)
		TheFormat = NH_LUMPED;
	else if( strcmp(Format,"Lined")==0  || strcmp(Format,"LINED")==0  || strcmp(Format,"lined")==0 )
		TheFormat = NH_LINED;
	else if( strcmp(Format,"NonLumped")==0  || strcmp(Format,"NONLUMPED")==0  || strcmp(Format,"Nonlumped")==0  || strcmp(Format,"nonlumped")==0)
		TheFormat = NH_NON_LUMPED;
	else {
		printf("\n\nSorry! Data Format %s not recognized!\n\nExiting to system...",Format);
	}
	
	/*  allocate memory to the TempLocTypes, and initialize to UNSET */
	TempLocTypes = (enum locus_type *)ECA_CALLOC((size_t)TempL,sizeof(enum locus_type));
	for(l=0;l<TempL;l++)  {
		TempLocTypes[l] = UNSET;
	}
	
	/*  get the first string on the next line.  If it is "LocusNames" then */
	/*  collect the locus names */
	fscanf(in,"%s",tempStr);
	if(strcmp(tempStr,"LocusNames") == 0)  {
		HasLocusNames = 1;
		TempLocNames = (char **)ECA_CALLOC((size_t)TempL,sizeof(char *));
		for(l=0;l<TempL;l++)  {
			TempLocNames[l]  = (char *)ECA_MALLOC(MAX_LOCUS_NAME_LENGTH * sizeof(char));
			fscanf(in,"%s",TempLocNames[l]);
		}
	}
	else {  /*  if it is not LocusNames, then we just ate the number of the first individual */
			/*  in the data set, and so we send that to the variable temp_int */
		temp_int = atoi(tempStr);
	}
	
	
	/*  allocate all necessary memory to the TempY array: */
	TempY = (int ***)ECA_CALLOC((size_t)TempM,sizeof(int **));
	for(i=0;i<TempM;i++) { double debugsilly;  /*  cycle over individuals, allocating necessary memory */
		TempY[i] = imatrix(TempL,2);
		/*debugsilly = gECA_mall_call_bytes; */
	}
	
	
	/* allocate memory to an array of arrays of ints that will ultimately be transferred over
	 to being the individual flags  */
	TempFlags = (int **)ECA_CALLOC((size_t)TempM, sizeof(int *));
	for(i=0;i<TempM;i++)  {
		TempFlags[i] = (int *)ECA_CALLOC((size_t)NUM_INDIV_FLAGS, sizeof(int));
	}
	
	/* allocate memory to an array of pointers to strings.  Leave them all NULL for now */
	TempIndNames = (char **)ECA_CALLOC((size_t)TempM, sizeof(char *));
	
	
	if(TheFormat == NH_LUMPED)  {
		ReadDataLumpedFormat(temp_int,HasLocusNames, TempY, 
					TempM, TempL, TempLocTypes, Divisor,  in, TempFlags, TempIndNames);
	}
	if(TheFormat == NH_NON_LUMPED) {
		ReadDataNonLumpedFormat(temp_int,HasLocusNames, TempY, 
					TempM, TempL, TempLocTypes,  in, TempFlags, TempIndNames);
	}
	
	
	
	
		
	/*  OK, now that this is all read into memory, we can start figuring out which of the loci */
	/*  are polymorphic and which locus indexes need shifting, etc.   */
	TempKl = ivector(0,TempL-1);
	TempIdx = imatrix(TempL,MAX_ALLELES+1);
	TempPolyL = 0;
	/* allocate some space for the allele names */
	TempAlleNames = (char ***)ECA_CALLOC(TempL,sizeof(char **));
	
	for(l=0;l<TempL;l++)  {
		TempKl[l] = 0;  /*  initialize to count the number of alleles */
		TempAlleNames[l] = (char **)ECA_CALLOC(MAX_ALLELES, sizeof(char *));
		for(i=0;i<TempM;i++)  {
			for(c=0;c<2;c++)  {
				if(TempY[i][l][c] >= 0 ) {	/*  only count alleles that are not, for example, "-1" (as those are missing data.		 */
					if(NotInArray(TempY[i][l][c],TempIdx[l],TempKl[l]))  {
						TempAlleNames[l][TempKl[l]] = (char *)ECA_CALLOC(MAX_ALLELE_NAME_LENGTH,sizeof(char));
						sprintf(TempAlleNames[l][TempKl[l]],"%d",TempY[i][l][c]);
						TempIdx[l][TempKl[l]++] = TempY[i][l][c];
					}
				}
			}
		}
		if(TempKl[l] > 1)
			TempPolyL++;
	}
	
	/*  we now have the resources to copy all that stuff over to the hyb_data */
	/*  First allocate the memory */
	D = (hyb_data *)ECA_MALLOC(sizeof(hyb_data));
	sprintf(D->DataFileName,"%s",FileName);
	
	/*  default the type of simulation to NewHybrids for now */
	D->TypeOfSim = NEW_HYBRIDS;
	D->StartingDispersionTheta = PRIORS_ONLY;
	D->StartingDispersionPi = PRIORS_ONLY;
	
	/* and some other defaults */
	D->PiTraceReport = 0;
	
	D->M = TempM;
	D->Kl = ivector(0,TempPolyL-1);
	D->Yobs = (ival ****)ECA_MALLOC(D->M * sizeof(ival ***));
	D->LocNames = (char **)ECA_MALLOC(TempPolyL * sizeof(char *));
	D->LocTypes = (enum locus_type *)ECA_MALLOC(TempPolyL * sizeof(enum locus_type));
	D->AlleleNames = (char ***)ECA_CALLOC(TempPolyL,sizeof(char **));
	for(l=0;l<TempPolyL;l++)  {
		D->LocNames[l] = (char *)ECA_MALLOC(MAX_LOCUS_NAME_LENGTH * sizeof(char));
	}
	for(i=0;i<D->M;i++)  {
		D->Yobs[i] = (ival ***)ECA_MALLOC(TempPolyL * sizeof(ival **));
		for(l=0;l<TempPolyL;l++)  {
			D->Yobs[i][l] = IvalVector(0,1,0,0,0);
		}
	}
	/*  Then copy over the polymorphic loci with the different indexes */
	for(i=0;i<D->M;i++)  {
		locus_counter = 0;
		for(l=0;l<TempL;l++)  {
			/*  pick out only those loci that have more than 1 allele */
			if(TempKl[l] > 1) {
				/* assign the string names of the alleles here */
				D->AlleleNames[locus_counter] = TempAlleNames[l];
				
				/*  set the Y values.  If the loci are Codom we jumble the */
				/*  names of the alleles however we wish, but for AFLP's the  */
				/*  O's must remain 0's and so with the 1's */
				for(c=0;c<2;c++)  {
					if(TempLocTypes[l] == CODOM)  {
						if(TempY[i][l][c] == -1)  {  /*  if it is missing data, set Yobs to missing too. */
							D->Yobs[i][locus_counter][c]->v = -1;
						}
						else {
							D->Yobs[i][locus_counter][c]->v = ReturnSubscript(
										TempIdx[l],TempKl[l],TempY[i][l][c]);
						}
					}
					else if(TempLocTypes[l] == AFLP) 
						D->Yobs[i][locus_counter][c]->v = TempY[i][l][c];
					else {
						printf("\n\nABORT!  We have a locus that is neither CODOM or AFLP\n\n");
						exit(1);
					}
						
				}
				locus_counter++;
			}
			/* if they are not polymorphic, then free the memory for the allele names
				actually, I am getting some weird memory errors.  So, I will comment this out
				and just plan on eating up a little extra memory.	*/
			/* else  {
				free(TempAlleNames[l]);
			} */
		}
	}
	
	/*  and now we should be able to do the allele numbers and locus numbers, types, and names, too */
	D->L = TempPolyL;
	D->AFLPsPresent = 0;
	locus_counter = 0;
	for(l=0;l<TempL;l++)  {
		if(TempKl[l] > 1) {
			/*  set the LocusNames and types */
			D->LocTypes[locus_counter] = TempLocTypes[l];
			if(D->LocTypes[locus_counter] == AFLP)
				D->AFLPsPresent = 1;
			if(HasLocusNames) 		/*  if they have names, then copy them */
				sprintf(D->LocNames[locus_counter],"%s",TempLocNames[l]);
			else  /*  if not, then give them a number designating which locus column they were from
					in the data set.  So, these may not go from 1 to the number of polymorphic loci */
				sprintf(D->LocNames[locus_counter],"%d",l+1);
			D->Kl[locus_counter] = TempKl[l];
			locus_counter++;
		}
	}
	
/*	
	out = fopen("ProcessedDataEcho.txt","w");
	fprintf(out,"%d  %d",D->M, D->L);
	for(i=0;i<D->M;i++)  {
		fprintf(out,"\n%d ",i+1);
		for(l=0;l<D->L;l++)  {
			fprintf(out,"  %d%d",D->Yobs[i][l][0]->v,D->Yobs[i][l][1]->v);
		}
	}
	fclose(out);
*/	
	
	/*  and finally, we must point the Flag variable of Data to TempFlags
		appropriate TempFlag array */
	D->IndFlags = TempFlags;
	/* then do the same with the names */
	D->IndNames = TempIndNames;
	
	/* and print out the allele name correspondences */
	out = fopen("aa-LociAndAlleles.txt","w");
	fprintf(out,"%d polymorphic loci from %d loci in data set",D->L,TempL);
	CYCLE_l(D)
		fprintf(out,"\n\nLocus %s",D->LocNames[l]);
		CYCLE_k(D) 
			fprintf(out,"\n\t%s",D->AlleleNames[l][k]);
		END1CYCLE
	END1CYCLE
	fclose(out);
	
	/*  now, free all the temporary memory */
	for(i=0;i<TempM;i++)  {
		free_imatrix(TempY[i],TempL);
	}
	free(TempY);
	free_ivector(TempKl,0,TempL-1);
	free_imatrix(TempIdx,TempL);
	free(TempLocTypes);
	if(HasLocusNames) {
		for(l=0;l<TempL;l++)  {
			free(TempLocNames[l]);
		}
		free(TempLocNames);
	}
	
	return(D);
}


/*  a function to read the genotype frequency class data from a file named FileName. */
/*  it allocates memory in D for this information */
void GetGtypFreqCats(hyb_data *D, char *FileName)
{
	int g, w[2];
	FILE *in, *out;
	
	/*  We can get all the gtyp freq data now too. */
	in = erdOpenFileOrRetry(FileName,"r");
	
	
	D->Gn = AllocIval(0,0,0);
	fscanf(in,"%d",&(D->Gn->v));
	/*  allocate memory to G and to CategoryNames */
	D->G = (double ***)ECA_MALLOC(D->Gn->v * sizeof(double **));
	D->CategoryNames = (char **)ECA_CALLOC((size_t)D->Gn->v,sizeof(char *));
	
	for(g=0;g<D->Gn->v;g++)  {
		D->G[g] = dmatrix(0,1,0,1);  /*  allocate memory */
		D->CategoryNames[g] = (char *)ECA_CALLOC(MAX_CAT_NAME_LENGTH, sizeof(char));
		/*  read the name of the category first */
		fscanf(in,"%s",D->CategoryNames[g] );
		/*  read the gtyp freq quantities */
		DO_W_CYCLE(w)
			fscanf(in,"%lf",&(D->G[g][w[0]][w[1]]) );
		END2CYCLE
	}
	fclose(in);
	
	
	/*  now print that to make sure it is ok */
	out = fopen("aa-EchoedGtypFreqCats.txt","w");
	for(g=0;g<D->Gn->v;g++)  {
		fprintf(out,"\n");
		DO_W_CYCLE(w)
			fprintf(out,"%.8f\t",D->G[g][w[0]][w[1]] );
		END2CYCLE
	}
	fclose(out);
	
}


/*  a function to copy over the GtypFreqCat information */
/*  from the command line options struct  */
void CopyGtypFreqCatsFromCL(hyb_data *D, cli_opts *CL)
{
	int g, w[2];
	FILE *out;
	
	D->Gn = AllocIval(0,0,0);
	D->Gn->v = CL->NumGProbs;
	
	/*  allocate memory to G and to CategoryNames */
	D->G = (double ***)ECA_MALLOC(D->Gn->v * sizeof(double **));
	D->CategoryNames = (char **)ECA_CALLOC((size_t)D->Gn->v,sizeof(char *));
	
	for(g=0;g<D->Gn->v;g++)  {
		D->G[g] = dmatrix(0,1,0,1);  /*  allocate memory */
		D->CategoryNames[g] = (char *)ECA_CALLOC(MAX_CAT_NAME_LENGTH, sizeof(char));
		
		sprintf(D->CategoryNames[g],"%s", CL->GProbNames[g]);
		
		/*  read the gtyp freq quantities */
		DO_W_CYCLE(w)
			D->G[g][w[0]][w[1]] = CL->GProbs[g][w[0]+w[1]];
			if(w[0]!=w[1]) {
				D->G[g][w[0]][w[1]] /= 2.0;  
			}
		END2CYCLE
	}
	
	
	/*  now print that to make sure it is ok */
	out = fopen("aa-EchoedGtypFreqCats.txt","w");
	for(g=0;g<D->Gn->v;g++)  {
		fprintf(out,"\n");
		DO_W_CYCLE(w)
			fprintf(out,"%.8f\t",D->G[g][w[0]][w[1]] );
		END2CYCLE
	}
	fclose(out);
	
}






/*  a function to read in the data when each individual's genotype is given */
/*  by a single number that must be split up. */
/*  temp_int is the first individual number if HasLocusNames == 0 */
/*  otherwise it just gets set to the individual number and then it */
/*  gets used throughout.  MISSING DATA IS DENOTED IN TEMPY BY -1 */
void ReadDataLumpedFormat(int temp_int, int HasLocusNames, int ***TempY, 
				int TempM, int TempL, enum locus_type *TempLocTypes, int Divisor, FILE *in, int **IndivFlags, char **IndNames)
{
	int i,l;
	char tempStr[100];
	FILE *out;
	int dummy;
	
	
	for(i=0;i<TempM;i++) {  /*  cycle over individuals */
		if(HasLocusNames==1 || i > 0) {
			fscanf(in," %s",tempStr);  /*  eat the first number (the fish index) */
			temp_int = atoi(tempStr);
		}
		if(temp_int != i+1)
			printf("\nEating the wrong fish number!!  Was expecting to find the number %d at the begninning of this line of data, but found %d instead.\n",i+1,temp_int);
			
		dummy = ReadIndivOptions(in, IndivFlags[i], i, &(IndNames[i]));
		
		/*  now cycle over loci and fill them in.  Get the locus thing first as a string */
		/*  then, after recording what type of locus it is, record it in Y. */
		for(l=0;l<TempL;l++)  {
			fscanf(in,"%s",tempStr);
			
			/*  the following function checks AND sets the locus_type */
			CheckLocType(tempStr,&TempLocTypes[l],i,l);
			if(TempLocTypes[l]==CODOM) {
				temp_int = atoi(tempStr);
				TempY[i][l][0] = temp_int / Divisor;
				TempY[i][l][1] = temp_int % Divisor;
				
				/*  call the locus missing data if either of the TempY's is zero */
				if(TempY[i][l][0] == 0 || TempY[i][l][1]==0)  {
					TempY[i][l][0] = -1;
					TempY[i][l][1] = -1;  
				}
			}
			else if(TempLocTypes[l]==AFLP) {
				if(strcmp(tempStr,"+")==0)  {
					TempY[i][l][0] = 1;
					TempY[i][l][1] = 1;
				}
				else if(strcmp(tempStr,"-")==0) {
					TempY[i][l][0] = 0;
					TempY[i][l][1] = 0;
				}
				else if(strcmp(tempStr,"0")==0) {
					TempY[i][l][0] = -1;
					TempY[i][l][1] = -1;
				}
				else 
					printf("\n\nError! Indiv %d, Locus %d.  AFLP but not + or -\n\n",i,l);
			}
			else if(TempLocTypes[l]==UNSET) {  /*  this will happen if the data are missing on
												   the first individual */
				if(strcmp(tempStr,"0")==0) {
					TempY[i][l][0] = -1;
					TempY[i][l][1] = -1;
				}
				else 
					printf("\n\nError! Indiv %d, Locus %d CheckLocType returned \nUNSET, but locus not missing!\n\n",i,l);
			}
			
		}
	}
	fclose(in);
	
	
	/*  print out the input to a file for comparison */
	out = fopen("EchoedGtypData.txt","w");
	for(i=0;i<TempM;i++) {  /*  cycle over individuals		 */
		fprintf(out,"\n%d",i+1);  /*  eat the first number (the fish index) */
	
		/*  now cycle over loci and fill them in */
		for(l=0;l<TempL;l++)  {
			fprintf(out,"  %d%d",TempY[i][l][0],TempY[i][l][1]);
			
		}
		
	}
	fclose(out);
	
}





/*  a function to read in the data when each individual's genotype is given */
/*  by a two numbers adjacent to each other on a line */
/*  temp_int is the first individual number if HasLocusNames == 0 */
/*  otherwise it just gets set to the individual number and then it */
/*  gets used throughout. */
void ReadDataNonLumpedFormat(int temp_int, int HasLocusNames, int ***TempY, 
				int TempM, int TempL, enum locus_type *TempLocTypes, FILE *in, int **IndivFlags, char **IndNames)
{
	int i,l, dummy;
	char tempStr[100];
	FILE *out;
	
	
	for(i=0;i<TempM;i++) {  /*  cycle over individuals */
		if(HasLocusNames==1 || i > 0) {
			fscanf(in," %s",tempStr);  /*  eat the first number (the fish index) */
			temp_int = atoi(tempStr);
		}
		if(temp_int != i+1)
			printf("\nEating the wrong fish number!! Was expecting to find the number %d at the begninning of this line of data, but found %d instead.\n",i+1,temp_int);
			
		dummy = ReadIndivOptions(in, IndivFlags[i], i, &(IndNames[i]));
		
		/*  now cycle over loci and fill them in.  Get the locus thing first as a string */
		/*  then, after recording what type of locus it is, record it in Y. */
		for(l=0;l<TempL;l++)  {
			fscanf(in,"%s",tempStr);
			
			/*  the following function checks AND sets the locus_type */
			/*  if it is codominant, it has to read the other allele as well. */
			/*  if it is AFLP it does not have to read the other allele and just */
			/*  records the info and proceeds */
			CheckLocType(tempStr,&TempLocTypes[l],i,l);
			if(TempLocTypes[l]==CODOM) {
				temp_int = atoi(tempStr);
				TempY[i][l][0] = temp_int;
				
				/*  then also grab the next allele */
				fscanf(in,"%s",tempStr);
				temp_int = atoi(tempStr);
				TempY[i][l][1] = temp_int;
				
				/*  then check to see if they were missing data.  Bark an error and exit it */
				/*  only one of the two was missing */
				if(TempY[i][l][0] <= 0 && TempY[i][l][1] <= 0)  {
					TempY[i][l][0] = -1;
					TempY[i][l][1] = -1;
				}
				else if(  (TempY[i][l][0] > 0 && TempY[i][l][1] <= 0) )  {
					printf("\n\nError!  Gene Copy 0 Non-Missing and Gene Copy 1 Missing");
					printf("\nAt Locus %d on Individual %d\n\n", l+1,  i+1);
					printf("Exiting to system...\n\n");
					exit(1);
				}
				
				else if(TempY[i][l][0] <= 0 && TempY[i][l][1] > 0) {
					printf("\n\nError!  Gene Copy 1 Non-Missing and Gene Copy 0 Missing");
					printf("\nAt Locus %d on Individual %d\n\n", l+1,  i+1);
					printf("Exiting to system...\n\n");
					exit(1);
				}
						   
			}
			else if(TempLocTypes[l]==AFLP) {
				if(strcmp(tempStr,"+")==0)  {
					TempY[i][l][0] = 1;
					TempY[i][l][1] = 1;
				}
				else if(strcmp(tempStr,"-")==0) {
					TempY[i][l][0] = 0;
					TempY[i][l][1] = 0;
				}
				else 
					printf("\n\nError! Indiv %d, Locus %d.  AFLP but not + or -\n\n",i,l);
			}
			
		}
	}
	fclose(in);
	
	
	/*  print out the input to a file for comparison */
	out = fopen("EchoedGtypData.txt","w");
	for(i=0;i<TempM;i++) {  /*  cycle over individuals		 */
		fprintf(in,"\n%d",i+1);  /*  print the first number (the fish index) */
	
		/*  now cycle over loci and fill them in */
		for(l=0;l<TempL;l++)  {
			fprintf(in,"   %d %d",TempY[i][l][0],TempY[i][l][1]);
			
		}
		
	}
	fclose(out);
	
}


/* this is a function that reviews all the data that was input for cline stuff */
void ProcessClineOptions(hyb_data *D, cline_struct *C)
{
	int i;
	int HasError=0;
	double minx=9999999999.0, maxx=-9999999999.0, lox, hix, span;
	
	
	if(C->NumLocales==0) { /* nothing to do in this case */
		return;
	}
	
	/* check to make sure we have locations for all locales and figure out the min and max distance */
	printf("Review of Cline Estimation Settings\n");
	for(i=0;i<=C->MaxLocaleIndex;i++) {
		if(C->HasLocale[i]>0) {
			printf("Locale Index %d has %d Indivs at location: ",i,C->NumInds[i]);
		}
		if(C->x[i] == LOCATION_UNSET) {
			printf("UNSET (This is an error! program will terminate)\n");
			HasError=1;
		}
		else {
			printf("%f\n",C->x[i]);
			if(C->x[i]<minx) {
				minx = C->x[i];
			}
			if(C->x[i] > maxx) {
				maxx = C->x[i];
			}
		}
	}
	
	if(HasError) {
		exit(1);
	}
	
	/* determine the points for drawing the clines */
	span = maxx - minx;
	lox = minx - .07 * span;
	hix = maxx + .07 * span;
	span = (hix-lox) / (C->NumClineXs - 1.0);  /* turn it into the increment */
	for(i=0,lox=lox;i<C->NumClineXs;i++,lox+=span) {
		C->ClineX[i] = lox;
	}
	
	/* now allocate memory to all the loci */
	C->Locs = (locline_struct **)ECA_CALLOC(D->L, sizeof(locline_struct *));
	for(i=0;i<D->L;i++) {
		C->Locs[i] = AllocLocline(C->MaxLocaleIndex,  C->NumClineXs,  10.0,  .005);
	}
	
}


locline_struct *AllocLocline(int MaxLocaleIndex,  double NumClineXs,  double asd,  double bsd)
{
	locline_struct *ret = (locline_struct *)ECA_MALLOC(sizeof(locline_struct));
	
	
	/* allocations */
	ret->AlleFreq = DvalVector(0,MaxLocaleIndex,  0,0,-1.0);
	ret->Y = IvalVector(0,MaxLocaleIndex,  0,0,-1.0);
	ret->N = IvalVector(0,MaxLocaleIndex,  0,0,-1.0);
	
	ret->alpha = AllocDval(0,0,-1.0);
	ret->beta = AllocDval(0,0,-1.0);
	
	/* initialize to starting values (should do this elsewhere, but for now...*/
/*	ret->alpha->v = 140.0;
	ret->beta->v  = .01; */
	
	ret->a_sd = AllocDval(0,0,0);
	ret->b_sd = AllocDval(0,0,0);
	ret->LogLike = AllocDval(0,0,0);
	
	ret->ClineY = DvalVector(0,NumClineXs,0.0,1.0,.01);
	
	/* variable initialization */
	ret->a_sd->v = asd;
	ret->b_sd->v = bsd;
	ret->LogLike->v = -999999999999999.999;
	
	
	return(ret);
	
}




















