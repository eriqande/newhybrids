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





/*
 *  RunWOG_func.c
 *  newhyb
 *
 *  Created by Eric C. Anderson on 5/24/07.
 *
 */



/* this just has a function that gets used to run all the WOG things */
int RunWithoutGraphics(hyb_chain *C, int DoAsBurnIn, int DoAsReps)
{
	int i,g;
	long EndSeed1, EndSeed2;
	time_t StartTime, CurrTime, LastTime, RealRepStartTime, DummyTime;
	double SecondsPerSweep, SecondsRemain;
	
	
	StartTime = time(NULL);
	LastTime = time(NULL);
	
	printf("\nStarting newhybs execution at %s with seeds %d %d\n\n",ctime(&StartTime),C->Seed1,C->Seed2);
	
	/* then start doing the burn-in reps */
	for(i=0;i<DoAsBurnIn;i++)  {

		SingleSweep(C);
		IncrementValues(C,1);

		if(i%300==0 && i > 0) {
			OutputHistograms(C);
			fprint_PofZ(C);
			fprint_AlleleAverages(C);
			fprint_PiAverages(C);
			fprint_UniPriPofZ(C);
			
		}
		
		CurrTime = time(NULL);
		if(difftime(CurrTime,LastTime) > 5.0)  {  /* if it has been 5 seconds since a report */
			LastTime = time(NULL);
			printf("\n\n\n\n\n\nNewHybrids started: %s", ctime(&StartTime));
			printf("\nData = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d",C->Dat->DataFileName, PriorType2str(C->Pri->PiPriType), 
				PriorType2str(C->Pri->ThetaPriType), C->Seed1, C->Seed2);
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
	
	RealRepStartTime = time(NULL);
	
	/* then do the actual data collection reps */
	for(i=0;i<DoAsReps;i++)  {
/*printf("\ni=%d",i); */
		SingleSweep(C);
/*printf("\n\tSwept");*/
		IncrementValues(C,1);
/*printf("\n\tInremented"); **/

		if(i%300==0 && i > 0) {
/* printf("\n\t********Starting Output***********\n"); */
			OutputHistograms(C);
			fprint_PofZ(C);
			fprint_AlleleAverages(C);
			fprint_PiAverages(C);
			fprint_UniPriPofZ(C);
		}
		
		CurrTime = time(NULL);
		
/* weird stuff to get it to work with microsoft */
/*DummyTime = CurrTime;
CurrTime = DummyTime;*/
/* printf("\n\tDifftime = %d",difftime(CurrTime,LastTime)); */
		if(difftime(CurrTime,LastTime) > 5.0)  {  /* if it has been 5 seconds since a report */
			LastTime = time(NULL);
			printf("\n\n\n\n\n\n\n\n\nNewHybrids started: %s", ctime(&StartTime));
			printf("\nData = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d",C->Dat->DataFileName, PriorType2str(C->Pri->PiPriType), 
				PriorType2str(C->Pri->ThetaPriType), C->Seed1, C->Seed2);
			printf("\nCurrently on Sweep %d of %d",i,DoAsReps);
			SecondsPerSweep = difftime(LastTime,RealRepStartTime)/(double)i;
			printf("\nSweeps per second = %f     Seconds per sweep = %f  ",1.0/SecondsPerSweep,SecondsPerSweep);
			SecondsRemain = (SecondsPerSweep * (double)(DoAsReps-i));;
			printf("\nShould finish Sweeps in %.0f seconds.  %.2f minutes. %.4f hours.  %.4f days.", 
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
	
	/* output the final bits of information */
	OutputHistograms(C);
	fprint_PofZ(C);
	fprint_AlleleAverages(C);
	
	/* then let the user know it has finished */
	CurrTime = time(NULL);
	printf("\n\nProgram New Hybrids completed at %s ",ctime(&CurrTime));
	printf("\nData = %s, PiPrior = %s, ThetaPrior = %s, Seeds = %d %d, NumBurnIn = %d, NumRepsAfterBurnIn = %d",
			C->Dat->DataFileName, PriorType2str(C->Pri->PiPriType), 
				PriorType2str(C->Pri->ThetaPriType), C->Seed1, C->Seed2,DoAsBurnIn,DoAsReps);
	getsd(&EndSeed1,&EndSeed2);
	printf("\nSeedsAtEnding = %ld %ld", EndSeed1, EndSeed2);
	
	printf("\n\nOutput is in the following files:");
	printf("\n\taa-LociAndAlleles.txt  --- A listing of loci and alleles read");
	printf("\n\taa-Pi.hist  ---  text file of a histogram representation for");
	printf("\n\t\tthe posterior distribution of Pi");
	printf("\n\taa-PofZ.txt  ---  for each individual in the sample, the posterior probability");
	printf("\n\t\tthat it falls into each of the genotype frequency categories");
	printf("\n\taa-Theta.hist  ---- a histogram representation of the posterior distributions");
	printf("\n\t\tof the allele frequencies in the different populations");
	printf("\n\n");
	
	return(0);
	
}
