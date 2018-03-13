#ifdef UN_EXTERN
	#define GLOB
#else 
	#define GLOB extern
#endif



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "MCTypesEtc.h"
#include "ECA_utilities.h"


#include "NewHybrids.h"



/*  converts a prior_type enum to a string */
const char *PriorType2str(enum prior_type loc)
{
	
	switch (loc) {
		case UNIFORM: return "UNIFORM";
		case JEFFREYS: return "JEFFREYS";
		case USER_DIRCH: return "USER_DIRCH";
		case FIXED_PRIOR: return "FIXED_PRIOR";
		default: return "!!Error-Unknown Prior Type!!";
	}
}


void fprintCL_Probs(FILE *out, cli_opts *CL)
{
	int i,j;
	
	fprintf(out,"COMM_LINE_OPTS: DataFilePath : %s\n",CL->DataFilePath);
	fprintf(out,"COMM_LINE_OPTS: GtypCatFilePath : %s\n",CL->GtypCatFilePath);
	fprintf(out,"COMM_LINE_OPTS: NumberOfGtypFreqCategories : %d\n",CL->NumGProbs);
	
	for(i=0;i<CL->NumGProbs;i++)  {
		fprintf(out,"COMM_LINE_OPTS: GtypFreqCategory : %d  :  %s  : ",i,CL->GProbNames[i]);
		for(j=0;j<3;j++)  {
			fprintf(out," %f ",CL->GProbs[i][j]);
		}
		fprintf(out,"\n");
	}
	
	fprintf(out,"COMM_LINE_OPTS: ThetaPriorType : %s\n",PriorType2str(CL->ThetaPriType));
	fprintf(out,"COMM_LINE_OPTS: AlleFreqPriorPath : %s\n",CL->AlleFreqPriorPath);

	fprintf(out,"COMM_LINE_OPTS: Seeds :  %ld  %ld  :  unset if they are equal to zero\n",CL->Seed1,CL->Seed2);

	fprintf(out,"COMM_LINE_OPTS: PiPrior :  %s",PriorType2str(CL->PiPriType));
	if(CL->PiPriType==FIXED_PRIOR) {
		fprintf(out," : ");
		for(i=0;i<CL->NumGProbs;i++)  {
			fprintf(out," %f ",CL->PiFixedValues[i]);
		}
	}
	if(CL->PiPriType==USER_DIRCH) {
		fprintf(out," : ");
		for(i=0;i<CL->NumGProbs;i++)  {
			fprintf(out," %f ",CL->PiUserDirchValues[i]);
		}
	}
	fprintf(out,"\n");
	
	fprintf(out,"COMM_LINE_OPTS: NumBurnIn :  %d\n",CL->NumBurnIn);
	fprintf(out,"COMM_LINE_OPTS: NumPostBurnIn :  %d\n",CL->NumPostBurnIn);
	
} 


void fprint_PofZ(hyb_chain *C)
{
	int i,g;
	FILE *out;
	char FILESTRING[10000];
	sprintf(FILESTRING,"%saa-PofZ.txt",gPWD);
	
	out = fopen(FILESTRING,"w");
	
	/*  print first row */
	fprintf(out,"%d_sweeps\tIndivName",C->Lat->Ind[0]->PofZ[0]->NumAved);
	CYCLE_g(C->Dat)
		fprintf(out,"\t%.3f/%.3f/%.3f/%.3f",C->Dat->G[g][0][0], 
					C->Dat->G[g][0][1], C->Dat->G[g][1][0], C->Dat->G[g][1][1]);
	END1CYCLE
	/*  then cycle over individuals, and within them cycle over the GtypFreq Cats and */
	/*  print out their posterior probs of PofZ */
	
	CYCLE_i(C->Dat)
		fprintf(out,"\n%d",i+1);
		if(C->Dat->IndNames[i]==NULL) {
			fprintf(out,"\tNoName");
		}
		else {
			fprintf(out,"\t%s",C->Dat->IndNames[i]);
		}
		CYCLE_g(C->Dat)
			fprintf(out,"\t%.5f",C->Lat->Ind[i]->PofZ[g]->Ave);
		END1CYCLE
	END1CYCLE
	
	fprintf(out,"\n");
	fclose(out);
}




void fprint_UniPriPofZ(hyb_chain *C)
{
	int i,g;
	FILE *out;
	char FILESTRING[10000];
	sprintf(FILESTRING,"%saa-ScaledLikelihood.txt",gPWD);
	
	out = fopen(FILESTRING,"w");
	
	/*  print first row */
	fprintf(out,"%d_sweeps",C->Lat->Ind[0]->UniPriPofZ[0]->NumAved);
	CYCLE_g(C->Dat)
		fprintf(out,"\t%.3f/%.3f/%.3f/%.3f",C->Dat->G[g][0][0], 
					C->Dat->G[g][0][1], C->Dat->G[g][1][0], C->Dat->G[g][1][1]);
	END1CYCLE
	/*  then cycle over individuals, and within them cycle over the GtypFreq Cats and */
	/*  print out their scaled likelihoods */
	
	CYCLE_i(C->Dat)
		fprintf(out,"\n%d",i+1);
		CYCLE_g(C->Dat)
			fprintf(out,"\t%.5f",C->Lat->Ind[i]->UniPriPofZ[g]->Ave);
		END1CYCLE
	END1CYCLE
	
	fprintf(out,"\n");
	fclose(out);
}


/*  This is a little function that prints out a list of the loci, and their respective alleles
	along with the current posterior mean of the allele frequency.  So, we get six columns:
	
	LocNumber  LocName   AlleleNumber   AlleleName    FreqInPop0  FreqInPop1 SDinPop0 SDinPop2
	
	and that should work just fine
	
*/
void fprint_AlleleAverages(hyb_chain *C)
{

	int l,k;
	FILE *out;
	hyb_data *D = C->Dat;
	double n0,n1;
	char FILESTRING[10000];
	
	sprintf(FILESTRING,"%saa-ThetaAverages.txt",gPWD);
	out = fopen(FILESTRING,"w");
	
	/*  print the number of completed sweeps */
	fprintf(out,"UponCompletionOf_%d_sweeps\n",C->Lat->Theta[0][0][0]->NumAved);
	
	/* print the header line */
	fprintf(out,"LocNumber LocName AlleleNumber AlleleName FreqInPop0 FreqInPop1 SDInPop0 SDinPop1\n");
	
	/* now cycle over the loci, and, within those, the alleles and print out */
	CYCLE_l(D) 
		CYCLE_k(D)
			n0 = C->Lat->Theta[0][l][k]->NumAved;
			n1 = C->Lat->Theta[1][l][k]->NumAved;
			
			fprintf(out,"%d\t%s\t%d\t%s\t%f\t%f\t%f\t%f\n", l,D->LocNames[l],
														k,D->AlleleNames[l][k],
														C->Lat->Theta[0][l][k]->Ave,
														C->Lat->Theta[1][l][k]->Ave,
														sqrt(C->Lat->Theta[0][l][k]->Var*n0),
														sqrt(C->Lat->Theta[1][l][k]->Var*n1));
		END1CYCLE
	END1CYCLE

	
	fclose(out);
}


void fprint_PiAverages(hyb_chain *C)
{
	int g;
	FILE *out;
	char FILESTRING[10000];
	
	sprintf(FILESTRING,"%saa-Pi.aves",gPWD);
	
	if( (out=fopen(FILESTRING,"w"))==NULL) {
		fprintf(stderr,"Error: Can't open file \"Pi.Aves\" to write output.\n");
		exit(1);
	}
	fprintf(out,"GtypFreqCategory\tPosteriorMeanPi\n");
	CYCLE_g(C->Lat) 
		fprintf(out,"%s\t%f\n",C->Dat->CategoryNames[g],C->Lat->Pi[g]->Ave);
	END1CYCLE 
	fclose(out);
	
}



/* this prints the current values of pi in the chain.  If PrintHeader is 1, then it just prints a header. */
void printPi_Trace(int Rep, hyb_chain *C, int PrintHeader) 
{
	int g;
	
	if(PrintHeader) {
		printf("PI_TRACE:Rep");
		CYCLE_g(C->Lat) 
		printf("\t%s",C->Dat->CategoryNames[g]);
		END1CYCLE
		printf("\n");
		
		return;
	}
	
	printf("PI_TRACE:%d",Rep);
	CYCLE_g(C->Lat) 
		printf("\t%f",C->Lat->Pi[g]->v);
	END1CYCLE
	printf("\n");
}

/* this prints the current values of pi in the chain.  If PrintHeader is 1, then it just prints a header. */
void printZ_Trace(int Rep, hyb_chain *C, int PrintHeader)
{
    int i,g;
    char tempStr[1000];
    
    if(PrintHeader) {
        printf("Z_TRACE:# Categories (first is 0): ");
        CYCLE_g(C->Lat)
        printf("%s ",C->Dat->CategoryNames[g]);
        END1CYCLE
        printf("\n");
        
        printf("Z_TRACE:IndIdx\tIndName\tRep\tZ\n");
        
        return;
    }
    
    CYCLE_i(C->Dat)
    printf("Z_TRACE:");
    if(C->Dat->IndNames == NULL || C->Dat->IndNames[i]==NULL) {
        sprintf(tempStr,"NA");
    } else {
        sprintf(tempStr,"%s", C->Dat->IndNames[i]);
    }
    printf("%d\t%s\t%d\t%d\n", i, tempStr, Rep, C->Lat->Ind[i]->Z->v);
    END1CYCLE
}











