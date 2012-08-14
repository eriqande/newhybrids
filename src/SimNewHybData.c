/*
	Collection of routines, including a function main() for simulating
	data under a model similar (or identical) to that used in NewHybrids.  

	Output is a datafile that can then be read by NewHybrids. 
	
	This uses some of the data structures from NewHybs (i.e. for reading and
	storing the genotype freq categories, etc.)
	
	
	Command line syntax:
		
	
	
		-i g [z] [s] num  --> make num individuals in category g having options z and s  (so 'z' and 
								's' are optional arguments  to the -i flag.
		-a nloc pa pb  --> make nloc AFLP loci with frequencies pa and pb  in spp. 0 and 1 respectively

		-c filename  -->  make codominant loci, the number of which and the allele frequencies at which
						  will be given in the file "filename"
						  
						  The file should have this format:
						  
						  nloci
						  num_alleles_1
						  p1 p2 p3 ... pk
						  q1 q2 q3 ... qk
						  num_alleles_2
						  p1 p2 p3 ... pk
						  q1 q2 q3 ... qk
						   
						  and so forth...
						  
						  
		-a and -i may be used multiple times.  -c can be used multiple times as well  That
		way you can mash together multiple files of allele freqs into a big data set. 
		
		Example:
		
		-a 10 .1 .9 -a 7 .95 .02 -a 20 .7 .2 -i 0 z s 15 -i 0 20 -i 1 z s 15  -c  codom_freqs.in

*/ 

#define UN_EXTERN

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "ECA_MemAlloc.h"
#include "MCTypesEtc.h"
#include "MathStatRand.h"
#include "ranlib.h"

/* #include "GFMCMC.h" 
#include "GFMCMC_DrawingUtilities.h" */

#include "NewHybrids.h"
#include "SimNewHybData.h" 

#ifdef __MWERKS__
	#ifdef macintosh
		#include <console.h>
	#endif
#endif


/* prototypes of functions in this file */
char SimAFLPLocus(double pa, double pb, double **g);
void GetCodomFreqs(char *FileName, double ***Freqs, int *AlleNums, int *NLoci);
void SimCodomLocus(double **Freqs, int K, double **g, char *String);
void usage(void);


/* some pre-defined constants  */
#define MAX_INDIV 10000
#define MAX_LOCI  500
#define MAX_CAT   50
#define MAX_CODOM_LOCI 500
#define MAX_CODLOC_FILES 20




int main(int argc,char** argv)
{
	int i,j,l,z,s,t,nn,mm;
	hyb_data *D = (hyb_data *)ECA_MALLOC(sizeof(hyb_data));
	int **N,n,g, TotalN=0;  /*  the number of the different types of individuals to simulate:
						N[g][0]  is the number of category g's with no z or s options
						N[g][1]  is the number of category g's having the z==g option
						N[g][2]  is the number of category g's having the z==g option
									AND which have the s option as well.  */
	int M=0,m;  /* the number of loci */
	double ***p;  /* the allele freqs at  the different loci.
						p[0][l][k] is the freq of the k-th allele at the l-th locus in spp. 0
						p[1][l][k] is the same in spp. 1. 
						for AFLP loci p[0][l][0] is the freq of the observable '+' allele
						in spp. 0.  and p[0][l][1] is not going to be defined */
	enum locus_type *LT;  /* the types of the loci CODOM, AFLP, etc */
	double q0,q1;  /* temporary, for getting the AFLP allele freqs off the command line */
	FILE *out, *truths;
	char boing;
	char **CodLocFileNames = (char **)ECA_CALLOC(MAX_CODLOC_FILES, sizeof(char *));
	int NumCodLocFiles = 0;
	int *NumCodAlleles = (int *)ECA_CALLOC(MAX_CODOM_LOCI, sizeof(int));
	int NumCodLoci = 0;
	double ***CodLociFreqs = (double ***)ECA_CALLOC(MAX_CODOM_LOCI, sizeof(double **));
	char tempStr[40];
	char GtypsFile[100];
	long seed1, seed2;
	FILE *sfile;  /* for the seeds */
	time_t CurrTime;
	
	
	#ifdef __MWERKS__
	#ifdef macintosh
		argc = ccommand(&argv);
	#endif
	#endif	

	/* set default for gtyp file */
	sprintf(GtypsFile,"TwoGensGtypFreq.txt");
	
	/* allocate memory */
	N = (int **)ECA_CALLOC((size_t)MAX_CAT, sizeof(int *));
	for(i=0;i<MAX_CAT;i++)  
		N[i] = (int *)ECA_CALLOC( (size_t)3, sizeof(int *));
	
	p = (double ***)ECA_CALLOC(2,sizeof(double **));
	for(j=0;j<2;j++)  {
		p[j] = (double **)ECA_CALLOC(MAX_LOCI,sizeof(int *));
		for(i=0;i<MAX_LOCI;i++)  {
			p[j][i] = (double *)ECA_CALLOC(1,sizeof(int));
		}
	}
	
	LT = (enum locus_type *)ECA_CALLOC(MAX_LOCI, sizeof(enum locus_type));
	
	
	/* process command line options */
	if(argc<=2) {
		usage();
	}
	
	i=1;
	while(i<argc){
		switch(argv[i][1]){
			case 'a':  /* include some AFLP loci  */
				m = atoi(argv[++i]);
				q0 = (double)atof(argv[++i]);
				q1 = (double)atof(argv[++i]);
				
				for(mm=M+m;M<mm;M++) {
					LT[M] = AFLP;
					p[0][M][0] = q0;
					p[1][M][0] = q1;
				}
				i++;
				break;
			case 'i':  /* specify some individuals */	 
				z=s=t=0; 
				g = atoi(argv[++i]);  /* first number following it is the hybrid category */
				i++;
				while(isdigit(argv[i][0])==0)  {  /* keep chewing through (and processing) */
											      /* options till you get to the number of indivs */
					if(argv[i][1] != '\0')
							printf("\nz or s option expected, but found a string longer than 1 in command line!\n\n");
					if(argv[i][0]=='z')
						z = 1;
					if(argv[i][0]=='s')
						s = 1;
					i++;
				}				 
				n = atoi(argv[i++]);   /* number of individuals of category g with options as specified */
										/* also incrementing argv to the next value here */
				
				/* Process the z and s options  */
				if(s==1 && z==0)  {
					printf("\n\nCommand line error!  s option requires z, too!!\n\n");
					printf("Exiting to system...\n\n");
					exit(1);
				}
				else {
					if(s==1)
						t = 2;
					if(z==1 && s==0)
						t = 1;
				}
				/* then put the result in the appropriate place in the N array */
				N[g][t] += n;
				TotalN += n;
				break;
			case 'c':
				CodLocFileNames[NumCodLocFiles] = (char *)ECA_CALLOC(80, sizeof(char));
				sprintf(CodLocFileNames[NumCodLocFiles],"%s",argv[++i]);
				GetCodomFreqs(CodLocFileNames[NumCodLocFiles++], CodLociFreqs, NumCodAlleles, &NumCodLoci);
				i++;
				break;
			case 'g':
				sprintf(GtypsFile,"%s",argv[++i]);
				i++;
				break;
			default:
				printf("\n\nUh-oh! Something was wrong on the command line!!\n\n");
				usage();    
		}
	}
	
	
	
	/* tell them that the program is launching and give the command line: */
	CurrTime = time(NULL);
	printf("Beginning execution of simdata_nh at %s ",ctime(&CurrTime));
	printf("\nUsing command line:  ");
	for(i=0;i<argc;i++) printf("%s ",argv[i]);
	printf("\n");
	
	GetGtypFreqCats(D, GtypsFile);
	
	
	
	/* deal with the random number seeds */
	if( (sfile = fopen("snhd_seeds","r"))==NULL) {
		seed1 = (long)time(NULL);
		seed2 = seed1/3;
		printf("\nNo file \"snhd_seeds\".  Using clock to set RNG: %d %d\n", seed1,seed2);
	}
	else {
		fscanf(sfile,"%d %d",&seed1,&seed2);
		printf("\nUsing contents of snhd_seeds to set RNG: %d %d\n", seed1,seed2);
		fclose(sfile);
	}
	setall(seed1,seed2);

	
	
	while( (out = fopen("SimData.out","w"))==NULL) {
		printf("\nCan't open the file.  Might be open in another application.\n\n");
		printf("Close it and hit return\n");
		scanf("%c", &boing);
	}
	
	while( (truths = fopen("SimDataGs.out","w"))==NULL) {
		printf("\nCan't open the SimDataGs.out.  Might be open in another application.\n\n");
		printf("Close it and hit return\n");
		scanf("%c", &boing);
	}
	

	
	
	/* print the preamble to the data file and the truth file*/
	fprintf(truths,"IndivNumber TrueGtyp Options\n");
	fprintf(out,"NumIndivs %d\nNumLoci %d\nDigits 3\nFormat Lumped\n\n",TotalN,M+NumCodLoci);
	for(g=0,nn=0;g<D->Gn->v;g++)  {
		for(t=0;t<3;t++)  {
			for(n=0;n<N[g][t];n++)  {
				fprintf(out,"\n%d\t",++nn);
				fprintf(truths,"%d\t%d\t",nn,g);
				if(t==0)
					fprintf(truths,"n\t\n");
				if(t==1) {
					fprintf(truths,"z\t\n");
					fprintf(out,"z%d ",g);
				}
				if(t==2) {
					fprintf(out,"z%ds ",g);
					fprintf(truths,"zs\t\n");
				}
				
				/* simulate the codominant loci first */
				for(l=0;l<NumCodLoci;l++) {
					SimCodomLocus(CodLociFreqs[l], NumCodAlleles[l], D->G[g], tempStr);
					fprintf(out, " %s",tempStr);
				}
				
				/* then put all the AFLP loci after them */
				for(l=0;l<M;l++)  {
					if(LT[l]==AFLP)
						fprintf(out," %c",SimAFLPLocus(p[0][l][0],p[1][l][0],D->G[g]) );
				}
				
			}
		}
	}
	
	fclose(out);
	fclose(truths);
	
	
	/* down here we will print the allele freqs out to a file called SimAlleFreqs.out
	   The format will be LocType: p1 p2 ... pK, where LocType is AFLP or CODOM */
	if( (out=fopen("SimAlleFreqs.out","w"))==NULL) {
		fprintf(stderr,"\nWARNING!  Can't open file SimAlleFreqs.out to write allele freqs to it\n",seed1,seed2);
	}
	else {
		/* report the codominant loci first */
		i=0;
		for(l=0;l<NumCodLoci;l++) {int k;
			fprintf(out,"CODOM: %d : ",++i);
			for(k=0;k<NumCodAlleles[l];k++)
				fprintf(out, "%f ",CodLociFreqs[l]);
			fprintf(out,"\n");
		}
		
		/* then put all the AFLP loci after them */
		for(l=0;l<M;l++)  {
			if(LT[l]==AFLP)
				fprintf(out,"AFLP: %d : %f %f\n",++i,p[0][l][0],p[1][l][0] );
		}
		fclose(out);
	}
	
	
	
	getsd(&seed1,&seed2);
	if( (sfile=fopen("snhd_seeds","w"))==NULL) {
		printf("\nWARNING!  Can't open file snhd_seeds to write seeds %d %d to it\n",seed1,seed2);
	}
	else {
		printf("\n\nFinished on random seeds: %d %d : printed to snhd_seeds",seed1,seed2);
		fprintf(sfile,"%d %d\n",seed1,seed2);
		fclose(sfile);
	}
	printf("\n\nAll Done\n\nSimulated Gtyps in File: SimData.out\n\nActual hybrid categories in file SimDataGs.out\n\n");
	
	return(0);
}

/*
	Given allele freqs for the presence of a band at an AFLP marker being 
	pa and pb in pop's 0 and 1, respectively, and given the hybrid category
	having the genotype frequencies as stored in g (as g[0][0], g[0][1],
	g[1][0], and g[1][1]) this function simulates the 
	observed state of the individual at this locus and returns either a '+' or  '-'.  
*/
char SimAFLPLocus(double pa, double pb, double **g)
{
	double PlusProb=0.0, rando;
	
	/* for the prob that it's 00  */
	PlusProb += g[0][0] * (1.0 - (1.0-pa)*(1.0-pa));
	
	/* for the prob that it's 11 */
	PlusProb += g[1][1] * (1.0 - (1.0-pb)*(1.0-pb)); 
	
	/* for the prob that it's 01 (note that g[0][1] and g[1][0] should be the same,
	but I compute these separately for generality  */
	PlusProb += g[0][1] * (1.0 - (1.0-pa)*(1.0-pb));
	
	/* finally we collect the prob for 10 */
	PlusProb += g[1][0] * (1.0 - (1.0-pb)*(1.0-pa));
	
	/* at this junction, PlusProb is the probability that the locus is a '+', so
	now all we have to do is draw the thing */
	rando = (double)ranf();
	
	if(rando<PlusProb)
		return('+');
		
	return('-');
	
}

/* for a locus with K alleles, with allele freqs in the two species given
   in Freqs, simulate a genotype given genotype freqs stored in g.
   Return it as a pointer to string. */
void SimCodomLocus(double **Freqs, int K, double **g, char *String) {
	int i,j, first, second;
	double rando = (double)ranf();
	double sum = 0.0;
	char Str1[10], Str2[10];
	
	/* first determine the gtyp freq category */
	for(i=0;i<2;i++)  {
		for(j=0;j<2;j++)  {
			sum += g[i][j];
			if (sum > rando) {
				first = i;
				second = j;
				i=j=2;  /* get out of the loop */
			}
			
		}
	}
	
	/* now store the allelic types in i and j.  Remember, they
	range from 1 to K, with 0 being reserved for missing values  */
	i = IntFromProbsRV(Freqs[first], 1, K);
	j = IntFromProbsRV(Freqs[second], 1, K);
	
	/* then form the strings, assuming NumDigits = 3 */
	if(i<10)
		sprintf(Str1,"00%d",i);
	else 
		sprintf(Str1,"0%d",i);
	if(j<10)
		sprintf(Str2,"00%d",j);
	else 
		sprintf(Str2,"0%d",j);
	
	sprintf(String, "%s%s", Str1, Str2);
}

void GetCodomFreqs(char *FileName, double ***Freqs, int *AlleNums, int *NLoci)
{
	FILE *in;
	int i,j,k,NumLoc, start;
	double sum;
	
	while( !(in = fopen(FileName,"r") ) ) {
		fprintf(stderr,"\nCan't open file %s. NLoci = %d\n\nExiting!\n\n",FileName,*NLoci);
		exit(1);
	}
	
	fscanf(in,"%d",&NumLoc);
	start = *NLoci;
	*NLoci += NumLoc;
	
	/* echo this stuff to stdin */
	printf("\n\n-----------------------\n\n");
	printf("Reading %d Loci from \"%s\" (from %d to %d)",NumLoc,FileName,start, *NLoci-1);
	
	/* now cycle over the loci */
	for(j=start;j<*NLoci;j++)  {
		fscanf(in,"%d",&AlleNums[j]);
		printf("\n%d", AlleNums[j]);
		Freqs[j] = (double **)ECA_CALLOC(2, sizeof(double *));
		for(i=0;i<2;i++)  {
			printf("\n\t");
			Freqs[j][i] = (double *)ECA_CALLOC((size_t)AlleNums[j],sizeof(double));
			sum = 0.0;
			for(k=0;k<AlleNums[j];k++)  {
				fscanf(in,"%lf", &Freqs[j][i][k]);
				printf("%.4f\t", Freqs[j][i][k]);
				sum += Freqs[j][i][k];
			}
			if(sum < .99999999 || sum > 1.0000001)  {
				fprintf(stderr, "\nWarning! Sum of allele freqs not 1.0.  Locus %d, in %s\n\n",
					j - start, FileName);
			}
		}
	}
	printf("\n\nDone with file \"%s\"\n_________________",FileName);
	
	
	fclose(in);
	
	
}


void usage(void)
{
	printf("\n\nThis is a program to simulate data under the newhybrids model.");
	printf("\nIt outputs data in newhybrids format.");
	printf("\nThe program looks for random number seeds in a file called \"snhd_seeds\".");
	printf("\nThe contents of that file should be two integers on a single line.");
	printf("\nIf that file is not found, a random seed is made from the current time.");
	printf("\nUpon completion, the program writes the current random number generator");
	printf("\nstate to the file snhd_seeds to be used on the next run.");
	printf("\n\nCommand line syntax:");
		printf("\n\t-g GtypFile   --> where GtypFile is the pathname to the genotype freq category file.");
		printf("\n\t-i g [z] [s] num  --> make num individuals in category g having options z and s  (so 'z' and");
		printf("\n\t						's' are optional arguments  to the -i flag.");
		printf("\n\t-a nloc pa pb  --> make 10 AFLP loci with frequencies pa and pb spp. 0 and 1 respectively");
		printf("\n\t-c filename  -->  make codominant loci, the number of which and the allele frequencies at which");
		printf("\n\t				  will be given in the file \"filename\"");
						  
		printf("\n\t				  The file should have this format: ");
		printf("\n");				  
		printf("\n\t				  nloci");
		printf("\n\t				  num_alleles_1");
		printf("\n\t				  p1 p2 p3 ... pk");
		printf("\n\t				  q1 q2 q3 ... qk");
		printf("\n\t				  num_alleles_2");
		printf("\n\t				  p1 p2 p3 ... pk");
		printf("\n\t				  q1 q2 q3 ... qk");
		printf("\n");				   
		printf("\n\t				  and so forth...");
						  
						  
		printf("\n\n\n-a and -i may be used multiple times.  -c can be used multiple times as well!  That");
		printf("\nway you can mash together multiple files of allele freqs into a big data set."); 
		
		printf("\n\nExample:");
		
		printf("\n\n-a 10 .1 .9 -a 7 .95 .02 -a 20 .7 .2 -i 0 z s 15 -i 0 20 -i 1 z s 15  -c  codom_freqs.in\n\n");
		
		exit(1);
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








