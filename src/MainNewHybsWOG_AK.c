#define UN_EXTERN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "ECA_MemAlloc.h"
#include "MCTypesEtc.h"
#include "ECA_utilities.h"
#include "ECA_ReadData.h"
#include "MathStatRand.h"
#include "ranlib.h"

#include "NewHybrids.h"


int main(void)
{	
	char File[100], DatFile[100], GtypFile[100];
	int temp;
	double temp_d;
	enum prior_type PiP = JEFFREYS, ThP = JEFFREYS;
	int DoAsBurnIn = 0, DoAsReps = 0,i,g;
	time_t StartTime, CurrTime, LastTime, RealRepStartTime, DummyTime;
	double SecondsPerSweep, SecondsRemain;
	int priorchoice;
	long Seed1, Seed2;
	int boing;
	
	hyb_data *D;
	hyb_prior *P;
	hyb_chain *C;
	
	printf("Welcome to the program currently named \"NewHybrids\" Version 1.1 beta");
	printf("\nPre-Released 2 OCTOBER 2002.");
	printf("\nwritten by Eric C. Anderson (eriq@u.washington.edu)");
	printf("\nCopyright (c) by The Regents of the University of California\n");
	printf("Please see user documentation for full software agreement.");
	
		
	printf("\n-----------------------\nEnter a \n\t0 to load \"TestDat.txt\"");
	printf("\n\t1 to load \"TestDatWithOptions.dat\"");
	printf("\n\t2 to load \"TestAFLP.dat\"");
	printf("\n\t3 to load \"TestAFLPWithOptions.dat\"");
	printf("\nas the Data File\n");
	printf("\nOr enter the name of the file you wish to use.");
	printf("\n\nRemember that file must be in the directory in which\nthe program resides.\n\n->");

	
	
	
	/*  read the DATA FILE name or a 0 */
	temp = erdGetNext(File,&temp_d,stdin);
	
	if(temp == 1 && temp_d == 0.0) {  /*  if input was an int equal to zero set File to "TestDat.txt" */
		sprintf(File,"TestDat.dat");
	}
	if(temp == 1 && temp_d == 1.0) { 
		sprintf(File,"TestDatWithOptions.dat");
	}
	if(temp == 1 && temp_d == 2.0) {  
		sprintf(File,"TestAFLP.dat");
	}
	if(temp == 1 && temp_d == 3.0) {  
		sprintf(File,"TestAFLPWithOptions.dat");
	}
	/*  get the DATA: */
	D = GetData(File);
	
	sprintf(DatFile,"%s",File);
	
	
	/*  now read in the gtyp freq categories */
	printf("\n\n-----------------\n");
	printf("\nNow enter a 0 to read the genotype frequency classes in the file");
	printf("\"TwoGensGtypFreq.txt\".  Or enter the name of the file holding your");
	printf("\nown definitions of genotype frequency classes.\n\n->");
	
	/*  read the GTYP FILE name or a 0 */
	temp = erdGetNext(File,&temp_d,stdin);
	
	if(temp == 1 && temp_d == 0.0) {  /*  if input was an int equal to zero set File to "TestDat.txt" */
		sprintf(File,"TwoGensGtypFreq.txt");
	}
	GetGtypFreqCats(D,File);
	sprintf(GtypFile,"%s",File);
	
	
	/*  now we can process the individual options a little bit  */
	ProcessIndivOptions(D);
	
	/* put alleles from indivs known to be in Purebred categories into the 
		"PriorizedAllelesPile" */
	PriorizeAllelesFromFixedZ(D);
	AddPriorizeAllelesFromFile(D);
	
	
	printf("\n\nGive me two small integers (>0) for random number seeds\n\n");
	scanf("%d%d",&Seed1,&Seed2);
	setall((long)Seed1,(long)Seed2);
	
	
	
	printf("\n\nEnter the choices for prior type for pi (mixing proportion):\n\t0 for Jeffreys\n\t1 for Uniform\n->");
	scanf("%d",&priorchoice);
	switch(priorchoice) {
		case(0): 
			PiP = JEFFREYS;
			break;
		case(1):
			PiP = UNIFORM;
			break;
		default:
			PiP = JEFFREYS;
			printf("\n\nInvalid choice.  Pi Prior set to Jeffreys by default.\n\n");
	}
	
	printf("\n\nEnter the choices for theta type for pi (mixing proportion):\n\t0 for Jeffreys\n\t1 for Uniform\n->");
	scanf("%d",&priorchoice);
	switch(priorchoice) {
		case(0): 
			ThP = JEFFREYS;
			break;
		case(1):
			ThP = UNIFORM;
			break;
		default:
			ThP = JEFFREYS;
			printf("\n\nInvalid choice.  Theta Prior set to Jeffreys by default.\n\n");
	}
		
	P = CreatePriors(D,PiP,ThP);
	C = CreateLatentChain(D,P);
	
	C->Seed1 = Seed1;
	C->Seed2 = Seed2;
		
	InitializeChain(C);
	
	printf("\n\nData all read and ready!!!\n\n");
	
	
	printf("\n\nEntering Graphics-Free version...\n\n");
	
	printf("\nPlease enter the number of sweeps for BurnIn (must be an integer) \n->");
	scanf("%d",&DoAsBurnIn);
	
	printf("\nPlease enter the number of sweeps to be done after BurnIn (must be an integer) \n->");
	scanf("%d",&DoAsReps);
	
	
	StartTime = time(NULL);
	LastTime = time(NULL);
	
	printf("\nStarting execution at %s\n\n",ctime(&StartTime));
	
	/* then start doing the burn-in reps */
	for(i=0;i<DoAsBurnIn;i++)  {

		SingleSweep(C);
		IncrementValues(C,1);

		if(i%300==0 && i > 0) {
			OutputHistograms(C);
			fprint_PofZ(C);
		}
		
		CurrTime = time(NULL);
		if(difftime(CurrTime,LastTime) > 5.0)  {  /* if it has been 5 seconds since a report */
			LastTime = time(NULL);
			printf("\n\n\n\n\n\nNewHybrids started: %s", ctime(&StartTime));
			printf("\nData = %s, GtypFile = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d",DatFile,GtypFile, PriorType2str(PiP), PriorType2str(ThP), Seed1, Seed2);
			printf("\nCurrently on BurnIn Sweep %d of %d",i,DoAsBurnIn);
			SecondsPerSweep = difftime(LastTime,StartTime)/(double)i;
			printf("\nSweeps per second = %f     Seconds per sweep = %f  ",1.0/SecondsPerSweep,SecondsPerSweep);
			SecondsRemain = (SecondsPerSweep * (double)(DoAsBurnIn-i));
			printf("\nShould finish BurnIn in %.0f seconds.  %.2f minutes.  %.4f hours.   %.4f days.", 
					SecondsRemain,SecondsRemain/60,SecondsRemain/3600,SecondsRemain/(3600*24)); 
					
			/* print out the current pi and the average pi */
			printf("\n\nCurrent Pi:\n");
			CYCLE_g(C->Dat)
				printf("  %.4f",C->Lat->Pi[g]->v);
			END1CYCLE
			
			printf("\nAverage Pi:\n");
			CYCLE_g(C->Dat)
				printf("  %.4f",C->Lat->Pi[g]->Ave);
			END1CYCLE
			printf("\n\n\n\n\n\n");
			fflush(stdout);
		}
	}
	
	printf("\n\n\n************************\n\nDone with %d Burn-In Sweeps\n\n\n*********************\n\n",DoAsBurnIn); 
	/* then reset the averages */
	ResetAllAveragesEtc(C);
	
	printf("\n\n\n\nAverages reset...\n\nStarting data accumulation reps...\n\n");
	
	
	/* then do the actual data collection reps */
	for(i=0;i<DoAsReps;i++)  {

		SingleSweep(C);

		IncrementValues(C,1);

		printf("%d ",i);

		if(i%300==0 && i > 0) {

			OutputHistograms(C);
			fprint_PofZ(C);
		
			printf("\n\n\n\n\n\n\n\n\nNewHybrids Proceeding:");
			printf("\nData = %s, GtypFile = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d",DatFile,GtypFile, PriorType2str(PiP), PriorType2str(ThP), Seed1, Seed2);
			printf("\nCurrently on Sweep %d of %d",i,DoAsReps);
			
			/* print out the current pi and the average pi */
			printf("\n\nCurrent Pi:\n");
			CYCLE_g(C->Dat)
				printf("  %.4f",C->Lat->Pi[g]->v);
			END1CYCLE
			
			printf("\nAverage Pi:\n");
			CYCLE_g(C->Dat)
				printf("  %.4f",C->Lat->Pi[g]->Ave);
			END1CYCLE
			printf("\n\n\n\n\n\n");
			fflush(stdout);
		}
		
	}
	
	/* output the final bits of information */
	OutputHistograms(C);
	fprint_PofZ(C);
	
	/* then let the user know it has finished */
	printf("\n\nProgram New Hybrids completed. ");
	printf("\nData = %s, GtypFile = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d",DatFile,GtypFile, PriorType2str(PiP), PriorType2str(ThP), Seed1, Seed2);

	printf("\n\nOutput is in the following files:");
	printf("\n\taa-LociAndAlleles.txt  --- A listing of loci and alleles read");
	printf("\n\taa-Pi.hist  ---  text file of a histogram representation for");
	printf("\n\t\tthe posterior distribution of Pi");
	printf("\n\taa-PofZ.txt  ---  for each individual in the sample, the posterior probability");
	printf("\n\t\tthat it falls into each of the genotype frequency categories");
	printf("\n\taa-Theta.hist  ---- a histogram representation of the posterior distributions");
	printf("\n\t\tof the allele frequencies in the different populations");
	printf("\n\n");
	
	printf("\n\n\nEnter any integer to exit\n->");
	scanf("%d",&boing);
	
	return(0);
}
