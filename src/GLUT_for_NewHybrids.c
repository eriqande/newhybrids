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


/* we can leave a bunch of this out if we don't want a GUI */
#ifndef COMPILE_NEW_HYB_WITH_NO_GUI
#include "GFMCMC.h"
#include "GFMCMC_DrawingUtilities.h"
#include "GLUT_for_NewHybrids.h"

/*  some menu identifiers */
#define ALLELE_RIGHT_CLICK 0
#define CATEGORY_RIGHT_CLICK 1

/*  a global variable to record the max number of alleles at any locus */
GLOB int gMaxAlleNum;
GLOB char **gAlleStrings;
/*  a global variable to record the sort order for Kullback-Leibler */
GLOB int *gLocusSortArray;

/* a global variable to hold the sliding trace struct for the complete data loglikelihood */
GLOB sliding_trace_struct *gCompleteDataLogLike;



/*
	
	THESE ARE OUR MANDATORY FUNCTION DEFINITIONS.  THE PROTOTYPES
	ARE FOUND IN GLUT_FOR_MCMC.h 
	
	Here we must give definitions for:
	
	void gfmUserDefd_ResetAllAverages(void)
	void gfmUserDefd_InitializeChain(void)
	void OneStepByGlobals(void)
	
*/

void gfmUserDefd_ResetAllAverages(void)
{
	ResetAllAveragesEtc(gC);
	
	/*  then initialize the sliding trace structs, too
		this is GLUT specific, so it is separate from the
		standard ResetAllAverages Function */
	gfduInitSlidingTraceToZero(gCompleteDataLogLike);
}

void gfmUserDefd_InitializeChain(void)
{
	long temp1, temp2;
	
	if(gUseSameSeeds == 1) {
		setall((long)gC->Seed1, (long)gC->Seed2);
	}
	else {
		getsd(&temp1, &temp2);
		gC->Seed1 = temp1;
		gC->Seed2 = temp2;
	}
	InitializeChain(gC);
	
	
}


void gfmUserDefd_OneStepByGlobals(void)
{
	
	SingleSweep(gC);
	IncrementValues(gC,1);
	/* update the sliding Trace */
	gfduIncrementSlidingTrace(gCompleteDataLogLike,gC->Lat->TheCompleteDataLogLike);
	
	if(gNumSweepsAfterBurnIn%300==0 && gNumSweepsAfterBurnIn > 0) {
		OutputHistograms(gC);
		fprint_PofZ(gC);
		fprint_PiAverages(gC);
		fprint_AlleleAverages(gC);
		fprint_UniPriPofZ(gC);
	}
		
}




/*  this is the function that defines what all the windows are */
#define gfnhCAT_PROBS_SETUP gsDRAW_FUNC(gfnhDrawCategoryProbs); gsNUM_COLOR_KEYS( &(gC->Dat->Gn->v)); gsCOLOR_KEYS(gC->Dat->CategoryNames); gsMIDDLE_CLICK_MENU(CATEGORY_RIGHT_CLICK); gsCOLOR_SCHEME(DEEP_BLUE);
#define gfnhCAT_UNIPRI_PROBS_SETUP gsDRAW_FUNC(gfnhDrawCategoryUniPriProbs); gsNUM_COLOR_KEYS( &(gC->Dat->Gn->v)); gsCOLOR_KEYS(gC->Dat->CategoryNames); gsMIDDLE_CLICK_MENU(CATEGORY_RIGHT_CLICK); gsCOLOR_SCHEME(DEEP_BLUE);

#define gfnhALLELE_FREQ_SETUP gsDRAW_FUNC(gfnhDrawAlleleFreqs); gsNUM_COLOR_KEYS( &(gMaxAlleNum)); gsCOLOR_KEYS(gAlleStrings); gsMIDDLE_CLICK_MENU(ALLELE_RIGHT_CLICK); gsCOLOR_SCHEME(FISHER_PRICE);

void gfmUserDefd_DefineWindows(void)
{

	int l;
	/*  before getting down with it, here I compute the max num of alleles at any locus */
	/*  and also make a string that says "Allele 1" etc. */
	gMaxAlleNum = 0;
	CYCLE_l(gC->Dat)
		if(gC->Dat->Kl[l] > gMaxAlleNum)
			gMaxAlleNum = gC->Dat->Kl[l];
	END1CYCLE
	
	/*  now make the strings */
	gAlleStrings = (char **)ECA_CALLOC((size_t)gMaxAlleNum,sizeof(char *));
	for(l=0;l<gMaxAlleNum;l++)  {
		gAlleStrings[l] = (char *)ECA_CALLOC(30,sizeof(char));
		sprintf(gAlleStrings[l],"Allele %d",l+1);
	}

	gsOUTPUT_FILE_PREFIX("NewHybrids");
	
		
	/*  go ahead and do the definitions: */
	gsCONSOLE("\"Info\" Window");
		gsDRAW_FUNC(gfnhDrawConsoleWindow);
		gsKEYBOARD_FUNC(gfnhConsoleKeys);
	

	gsNEW_WINDOW(1,"Current Category Probabilities");
		gfnhCAT_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
			
	
	gsNEW_WINDOW(2,"Average Category Probabilities");
		gfnhCAT_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	
	gsNEW_WINDOW(3,"Current Allele Freqs");
		gfnhALLELE_FREQ_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
		
	gsNEW_WINDOW(4,"Average Allele Freqs");
		gfnhALLELE_FREQ_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	gsNEW_WINDOW(5,"Allele Frequency Histogram");
		gsDRAW_FUNC(gfnhDrawAlleleFreqsHist);
		gsCOLOR_KEYS(gAlleStrings);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsAXES_POS_UNIT_X;
		/*  note that NumColorKeys gets set depending on the Displayed Locus from */
		/*  within the function gfnhDrawAlleleFreqsHist() */
		
	gsNEW_WINDOW(6,"Category Probs Histogram");
		gsDRAW_FUNC(gfnhDrawCatProbsHist);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsAXES_POS_UNIT_X;
		
	gsNEW_WINDOW(7,"Observed Data");
		gsDRAW_FUNC(gfnhDrawObservedData);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		
		gsPADDING(-.17,1.17,-.17,1.17);
		
	
	gsNEW_WINDOW(11,"Current Origins Of Alleles");
		gsDRAW_FUNC(gfnhDrawWs);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		
		gsPADDING(-.17,1.17,-.17,1.17);
		
		
	gsNEW_WINDOW(12,"Average Origins Of Alleles");
		gsDRAW_FUNC(gfnhDrawW_Aves);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		gsPADDING(-.17,1.17,-.17,1.17);
		
		
		
	gsNEW_WINDOW(8,"Complete Data LogL Trace");
		gsDRAW_FUNC(gfnhDrawCompleteDataLogLTrace);
		gsCOLOR_SCHEME(FISHER_PRICE);
		
	gsNEW_WINDOW(10,"Current Kullback Leibler Div By Locus")
		gsDRAW_FUNC(gfnhDrawKullbLeib);
		gsXLO(-.15 * 12.0f);
		gsXHI(1.15 * 12.0f);
		gsYLO(-.15f * gC->Dat->L);
		gsYHI(1.15f * gC->Dat->L);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
		gsAXES_POS_QUADRANT
		gsNO_Y_AXIS
		
	gsNEW_WINDOW(9,"Average Kullback Leibler Div By Locus")
		gsDRAW_FUNC(gfnhDrawKullbLeib);
		gsXLO(-.15 * 12.0f);
		gsXHI(1.15 * 12.0f);
		gsYLO(-.15f * (GLfloat)gC->Dat->L);
		gsYHI(1.15f * (GLfloat)gC->Dat->L);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		gsAXES_POS_QUADRANT
		gsNO_Y_AXIS
		
	gsNEW_WINDOW(13, "Pi Values")
		gsDRAW_FUNC(gfnhDrawPiValues);
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		gsYLO(-2)
		gsYHI(2)
		gsXLO(-.6)
		gsXHI(1.15)
	
	gsNEW_WINDOW(14, "Pi Histogram")
		gsDRAW_FUNC(gfnhDrawPiHists);
		gsAXES_POS_QUADRANT;
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		
	gsNEW_WINDOW(15,"Current Scaled Likelihood");
		gfnhCAT_UNIPRI_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
			
	
	gsNEW_WINDOW(16,"Average Scaled Likelihoods");
		gfnhCAT_UNIPRI_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	gsNEW_WINDOW(20, "Allele Freq Clines");
		gsDRAW_FUNC(gfnhDrawAlleFreqClines);
		gsXLO( gClines->ClineX[0] - .1 * (gClines->ClineX[gClines->NumClineXs-1] - gClines->ClineX[0])  ); 
		gsXHI( gClines->ClineX[gClines->NumClineXs-1] );
		gsYLO(-.20);
		gsYHI(1.20);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsY_AXIS_CROSSES(0.0f);
		gsX_AXIS_CROSSES(0.0f);
		

}


void gfmUserDefd_DefineMenus(void)
{
	gUserDefdMenus[ALLELE_RIGHT_CLICK] = glutCreateMenu(gfnhOpenAlleleRightClickWindow);
	glutAddMenuEntry("Open Histogram View For Selected Locus", 5);


	gUserDefdMenus[CATEGORY_RIGHT_CLICK] = glutCreateMenu(gfnhOpenCategoryRightClickWindow);
	glutAddMenuEntry("Open Histogram View For Selected Individual", 0);
}





/*  this is a place that a user can do some last minute things like  */
/*  make some new menus and the like.  Here I don't have anything to do. */
void gfmUserDefd_LastWords(void)
{
	return;
}

void gfnhDrawConsoleWindow(void)
{
	char temp[250];
	
	glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Text);
	gfmRenderBitmapString(-.1f,1.0f,GLUT_BITMAP_HELVETICA_12, "Program: NewHybrids Version 1.1");
	gfmRenderBitmapString(-.1f,.9f,GLUT_BITMAP_HELVETICA_12, "Author: E.C. Anderson");
	
	sprintf(temp,"Data: \"%s\"", gC->Dat->DataFileName);
	gfmRenderBitmapString(-.1f,.8f,GLUT_BITMAP_HELVETICA_12, temp);
	
	sprintf(temp,"%d Sweeps",gNumSweeps);
	gfmRenderBitmapString(-.1f,.7f,GLUT_BITMAP_HELVETICA_12, temp);
	
	sprintf(temp,"%d After Burn In",gNumSweepsAfterBurnIn);
	gfmRenderBitmapString(-.1f,.6f,GLUT_BITMAP_HELVETICA_12, temp);
	
	sprintf(temp,"Seeds: %d   %d",gC->Seed1, gC->Seed2);
	gfmRenderBitmapString(-.1f,.5f,GLUT_BITMAP_HELVETICA_12, temp);
	
	if(gGo) 
		sprintf(temp,"RUNNING!");
	else 
		sprintf(temp,"STOPPED!");
		
	gfmRenderBitmapString(-.1f,.4f,GLUT_BITMAP_HELVETICA_12, temp);
	
	if(gC->Pri->ThetaPriType==JEFFREYS)
		sprintf(temp,"JEFFREYS");
	else if(gC->Pri->ThetaPriType==UNIFORM)
		sprintf(temp,"UNIFORM");
	else
		sprintf(temp,"UNKNOWN");
	sprintf(temp,"%s Prior for Theta (hit \"t\" to change)",temp);
	gfmRenderBitmapString(-.1f,.3f,GLUT_BITMAP_HELVETICA_10, temp);
	
	if(gC->Pri->PiPriType==JEFFREYS)
		sprintf(temp,"JEFFREYS");
	else if(gC->Pri->PiPriType==UNIFORM)
		sprintf(temp,"UNIFORM");
	else
		sprintf(temp,"UNKNOWN");
	sprintf(temp,"%s Prior for Pi (hit \"p\" to change)",temp);
	gfmRenderBitmapString(-.1f,.2f,GLUT_BITMAP_HELVETICA_10, temp);
	
	
}

/* handle the keyboard input to the console window, allowing 't'
   to toggle between Jeffreys and Uniform priors on Theta and
   'p' to do the same for Pi.  This overrides whatever 'p' and 't'
   might be in the Normal Key Processing Function.  */
void gfnhConsoleKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case('t'):
			SetThetaPriors(gC->Dat, gC->Pri, (enum prior_type)((int)(gC->Pri->ThetaPriType + 1)%2) );
			return;
		case('p'):
			SetPiPriors(gC->Dat, gC->Pri, (enum prior_type)((int)(gC->Pri->PiPriType + 1)%2) );
			return;
		default:
			gfmProcessNormalKeys(key, x, y);
	}
}

/*  draws histograms of each of the Pi's, uses different colors for JA_PURE and JA_ADMIXED */
void gfnhDrawPiHists(void)
{
	int g;
	gsCS;
	
	CYCLE_g(gC->Lat)
		glColor3fv(CS->Series[g]);
		gfmduDrawOneHistAsCurveWith2Bars(gC->Lat->Pi[g],1,GFMDU_NO_STIPPLE,0);
	END1CYCLE
	
	/*  then draw the axes */
	gfmDrawXYAxes();
	
}

/* draws the current and the average pi values as little bars  */
void gfnhDrawPiValues(void)
{
	int g;
	float avex=0.f,cx=0.f;
	gsCS
	gsWoH
	
	
	/* draw a series of overlapping rectangles going from 0 to 1 in x and from 
	1 to .5 and -1 to -.5 in y */
	CYCLE_g(gC->Lat)
		glColor3fv(CS->Series[g]);
		glRectf(cx,.5f,1.f,1.f);
		glRectf(avex,-.5f,1.f,-1.f);
	
		cx  += (float)gC->Lat->Pi[g]->v;
		avex +=(float)gC->Lat->Pi[g]->Ave;
	END1CYCLE
	
	/* then draw the text "Aves" and "Current"  */
	glColor3fv(CS->Text);
	gfmStrokeString("Current", .41f, -.05f, .5f, 4, 0.0f, WoH);
	gfmStrokeString("Average", .41f, -.05f, -1.0f, 4, 0.0f, WoH);
	
}

void gfnhDrawKullbLeib(void)
{
	int l;
	GLfloat x=0.0f,y,h=.7f,w,max=0.0f;
	char S[100];
	FILE *out;
	gsSettings gsCS gsWoH;
	
	
	CYCLE_l(gC->Dat)
		/* set the height and the width */
		y = (GLfloat)(gC->Dat->L - l);
		
		if(Settings->GFNH_CURRENT_OR_AVES == GFNH_AVES)
			w = (GLfloat)gC->Lat->Locus_KB[l]->Ave;
		else if(Settings->GFNH_CURRENT_OR_AVES == GFNH_CURRENT)
			w = (GLfloat)gC->Lat->Locus_KB[l]->v;
		else {
			out = fopen(gERROR_FILE,"a");
			fprintf(out,"\nGFNH_CURRENT_OR_AVES set to neither GFNH_AVES or");
			fprintf(out,"\nGFNH_CURRENT in gfnhDrawKullbLeib!! \n\n");
			fclose(out);
		}
		if(w>max)
			max = w;
		
		/*  set the string of the locus and print it */
		sprintf(S,"%s ",gC->Dat->LocNames[l]);
		glColor3fv(CS->Text);
		gfmStrokeString(S, h*.85f ,x, y-h, 4, 0.0f, WoH);
		
		/*  then print the rectangle in the 3rd color of the series */
		glColor3fv(CS->Series[5]);
		glRectf(x,y,x+w,y-h);
	END1CYCLE
	
	gsHANDLE_EXTREMA(0,max,0,gC->Dat->L);
	
	gfmDrawXYAxes();
}


void gfnhDrawCompleteDataLogLTrace(void)
{		
	glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Series[0]);
	
	gfduDrawSlidingTrace(gCompleteDataLogLike);
	
	gfmDrawXYAxes();

}



/*  draw the data as two rows for each individual.  Each gene copy fits into */
/*  a unit square (which may get stretched) and will be color coded.  Blank */
/*  squares of the axis color denote missing data. */
/*  Horizontal lines between the individuals are toggled using the Xaxis DrawIt variable */
void gfnhDrawObservedData(void)
{
	int i,l,c;
	GLfloat xs = .8f, ys = .8f;
	GLfloat x,y, yy;
	ColorScheme3f *CS = gWindowsSettings[glutGetWindow()]->ColorScheme;
	XY_axis_struct *XA = gWindowsSettings[glutGetWindow()]->Xaxis;
	clipping_struct *C = gWindowsSettings[glutGetWindow()]->Clips;
	GLfloat WoH = gfmWoH();
	char temp[100];
	GLfloat texth;
	
	/*  lay down the locus names as rotated text at the top of each column */
	/*  set the top height of the bars */
	y = 3.0f * (GLfloat)gC->Dat->M;
	texth = xs * (C->yhi - C->ylo) / (C->xhi - C->xlo) ;  /*  for textheight */
	CYCLE_l(gC->Dat)  /*  cycle over the loci */
		x = (GLfloat)l;  /*  set the "cursor" to the left edge of each locus rectangle */
		
		sprintf(temp," %s",gC->Dat->LocNames[l]);
		glColor3fv(CS->Text);
		gfmStrokeString(temp,texth,x+.5f*xs,y,1,-45.0f,WoH);
	END1CYCLE
	
	
	CYCLE_i(gC->Dat)
		y = 3.0f * (GLfloat)gC->Dat->M - (3.0f * (GLfloat)i);
		
		/*  draw text telling the individual's number */
		glColor3fv(CS->Text);
		sprintf(temp,"%d",i+1);
		
		gfmStrokeString(temp, 1.9f ,-.1f, y-1.8f, 4, 0.0f,WoH);
		CYCLE_l(gC->Dat)
			x = (GLfloat)l;
			yy = y;
			if(gC->Dat->LocTypes[l] == CODOM) { /*  draw alleles if they are codominant */
				for(c=0;c<2;c++)  {
					if(gC->Dat->Yobs[i][l][c]->v >= 0)  {  /*  if it is not missing data draw the square */
						glColor3fv(CS->Series[gC->Dat->Yobs[i][l][c]->v % CS->N]);  /*  set the color */
						glRectf(x,yy,x+xs,yy-ys);
					}
					else {  /*  if it is missing, then draw an empty box in axis color */
						glColor3fv(CS->Axes);  /*  set the color */
						glBegin(GL_LINE_LOOP);
							glVertex2f(x,yy);
							glVertex2f(x+xs,yy);
							glVertex2f(x+xs,yy-ys);
							glVertex2f(x,yy-ys);
						glEnd();
					}
					/*  move the second allele down one unit: */
					yy -= 1.0f;
				}
			}
			/*  now, if the locus is an AFLP locus, we draw a plus or a zero */
			if(gC->Dat->LocTypes[l] == AFLP) {
				if(gC->Dat->Yobs[i][l][0]->v == 1)  {
					sprintf(temp,"+");
					glColor3fv(CS->Series[0]);
					gfmStrokeString(temp,3.0f,x + (.5f * xs),y-2.4f + .3f * ys, 2,0.0f,WoH);
				}
				else if(gC->Dat->Yobs[i][l][0]->v == 0) {
					sprintf(temp,"-");
					glColor3fv(CS->Series[1]);
					gfmStrokeString(temp,2.2f,x + (.5f * xs),y-2.2f + .3f * ys, 2,0.0f,WoH);
				}
				else if(gC->Dat->Yobs[i][l][0]->v == -1) {  /* draw an axis-colored, solid rectangle if it is missing */
					glColor3fv(CS->Axes);  /*  set the color */
					glRectf(x,yy,x+xs,yy-1.0-ys);
				}
					
				
				
			}
		END1CYCLE
		
		/*  draw a separating line  */
		if(XA->DrawIt==1)  {
			glColor3fv(CS->Axes);
			glBegin(GL_LINES);
				glVertex2f(0.0f,y-2.5f);
				glVertex2f((GLfloat)gC->Dat->L-(1.0f-xs),y-2.5f);
			glEnd();
		}
			
	END1CYCLE
	
	gsHANDLE_EXTREMA(0.0f,(GLfloat)gC->Dat->L,0.0f, 3.0f * (GLfloat)gC->Dat->M);
	
	
}



/*  draw the W's as two rows for each individual.  Each gene copy's W fits into */
/*  a unit square (which may get stretched) and will be color coded according to which 
/*  population it comes from.  The ones with missing data should get a little bracket */
/*  in the  axis around them.  The Ws will be sorted so that gene copies from population 0 always
	appear on top.  */
/*  Horizontal lines between the individuals are toggled using the Xaxis DrawIt variable */
void gfnhDrawWs(void)
{
	int i,l,c;
	GLint a[2],tempa;
	GLfloat xs = .8f, ys = .8f;
	GLfloat extra = .05;  /* for drawing lines around missing data boxes */
	GLfloat x,y, yy;
	ColorScheme3f *CS = gWindowsSettings[glutGetWindow()]->ColorScheme;
	XY_axis_struct *XA = gWindowsSettings[glutGetWindow()]->Xaxis;
	clipping_struct *C = gWindowsSettings[glutGetWindow()]->Clips;
	GLfloat WoH = gfmWoH();
	char temp[100];
	GLfloat texth;
	
	/*  lay down the locus names as rotated text at the top of each column */
	/*  set the top height of the bars */
	y = 3.0f * (GLfloat)gC->Dat->M;
	texth = xs * (C->yhi - C->ylo) / (C->xhi - C->xlo) ;  /*  for textheight */
	CYCLE_l(gC->Dat)  /*  cycle over the loci */
		x = (GLfloat)l;  /*  set the "cursor" to the left edge of each locus rectangle */
		
		sprintf(temp," %s",gC->Dat->LocNames[l]);
		glColor3fv(CS->Text);
		gfmStrokeString(temp,texth,x+.5f*xs,y,1,-45.0f,WoH);
	END1CYCLE
	
	
	CYCLE_i(gC->Dat)
		y = 3.0f * (GLfloat)gC->Dat->M - (3.0f * (GLfloat)i);
		
		/*  draw text telling the individual's number */
		glColor3fv(CS->Text);
		sprintf(temp,"%d",i+1);
		
		gfmStrokeString(temp, 1.9f ,-.1f, y-1.8f, 4, 0.0f,WoH);
		CYCLE_l(gC->Dat)
			x = (GLfloat)l;
			yy = y;
			if(1) { /*  draw alleles if they are codominant OR if they are AFLP */
				/* first order the Ws */
				a[0] = gC->Lat->Ind[i]->W[l][0]->v;
				a[1] = gC->Lat->Ind[i]->W[l][1]->v;
				/*if(a[1]==0) {
					tempa=a[0];
					a[0] = a[1];
					a[1]=tempa;
				}*/
				for(c=0;c<2;c++)  {
					if(gC->Dat->Yobs[i][l][c]->v >= 0)  {  /*  if it is missing or not data draw the square */
						glColor3fv(CS->Series[a[c] % CS->N]);  /*  set the color */
						glRectf(x,yy,x+xs,yy-ys);
					}
					if(gC->Dat->Yobs[i][l][c]->v == -1) {  /*  if it is missing, then draw an empty box in axis color, slightly larger than the allele box */
						glColor3fv(CS->Axes);  /*  set the color */
						glBegin(GL_LINE_LOOP);
							glVertex2f(x,yy);
							glVertex2f(x+xs,yy);
							glVertex2f(x+xs,yy-ys);
							glVertex2f(x,yy-ys);
						glEnd();
					}
					/*  move the second allele down one unit: */
					yy -= 1.0f;
				}
			}
			
		END1CYCLE
		
		/*  draw a separating line  */
		if(XA->DrawIt==1)  {
			glColor3fv(CS->Axes);
			glBegin(GL_LINES);
				glVertex2f(0.0f,y-2.5f);
				glVertex2f((GLfloat)gC->Dat->L-(1.0f-xs),y-2.5f);
			glEnd();
		}
			
	END1CYCLE
	
	gsHANDLE_EXTREMA(0.0f,(GLfloat)gC->Dat->L,0.0f, 3.0f * (GLfloat)gC->Dat->M);
	
	
}



/*  draw the W's as two rows for each individual.  Each gene copy's W fits into */
/*  a unit square (which may get stretched) and will be color coded according to which 
/*  population it comes from.  The ones with missing data should get a little bracket */
/*  in the  axis around them.  The Ws will be sorted so that gene copies from population 0 always
	appear on top.  */
/*  Horizontal lines between the individuals are toggled using the Xaxis DrawIt variable */
void gfnhDrawW_Aves(void)
{
	int i,l,c;
	GLfloat xs = .8f, ys = .8f;
	GLfloat x,y, yy;
	ColorScheme3f *CS = gWindowsSettings[glutGetWindow()]->ColorScheme;
	XY_axis_struct *XA = gWindowsSettings[glutGetWindow()]->Xaxis;
	clipping_struct *C = gWindowsSettings[glutGetWindow()]->Clips;
	GLfloat WoH = gfmWoH();
	char temp[100];
	GLfloat texth;
	
	/*  lay down the locus names as rotated text at the top of each column */
	/*  set the top height of the bars */
	y = 3.0f * (GLfloat)gC->Dat->M;
	texth = xs * (C->yhi - C->ylo) / (C->xhi - C->xlo) ;  /*  for textheight */
	CYCLE_l(gC->Dat)  /*  cycle over the loci */
		x = (GLfloat)l;  /*  set the "cursor" to the left edge of each locus rectangle */
		
		sprintf(temp," %s",gC->Dat->LocNames[l]);
		glColor3fv(CS->Text);
		gfmStrokeString(temp,texth,x+.5f*xs,y,1,-45.0f,WoH);
	END1CYCLE
	
	
	CYCLE_i(gC->Dat)
		y = 3.0f * (GLfloat)gC->Dat->M - (3.0f * (GLfloat)i);
		
		/*  draw text telling the individual's number */
		glColor3fv(CS->Text);
		sprintf(temp,"%d",i+1);
		
		gfmStrokeString(temp, 1.9f ,-.1f, y-1.8f, 4, 0.0f,WoH);
		CYCLE_l(gC->Dat)
			x = (GLfloat)l;
			yy = y;
			if(1) { /*  draw alleles if they are codominant OR if they are AFLP */
				
				for(c=0;c<2;c++)  {
					if(gC->Dat->Yobs[i][l][c]->v >= 0)  {  /*  if it is not missing draw the square */
						glColor3fv(CS->Series[0]);  /*  set the color for the whole underlying rectangle */
						glRectf(x,yy,x+xs,yy-ys);
						glColor3fv(CS->Series[1]);  /* then set the color for the Species=1 part. */
						glRectf(x+xs*(1.0-gC->Lat->Ind[i]->W[l][c]->Ave),yy,x+xs,yy-ys);

					}
					if(gC->Dat->Yobs[i][l][c]->v == -1) {  /*  if it is missing, then draw an empty box in axis color, slightly larger than the allele box */
						glColor3fv(CS->Axes);  /*  set the color */
						glBegin(GL_LINE_LOOP);
							glVertex2f(x,yy);
							glVertex2f(x+xs,yy);
							glVertex2f(x+xs,yy-ys);
							glVertex2f(x,yy-ys);
						glEnd();
					}
					/*  move the second allele down one unit: */
					yy -= 1.0f;
				}
			}
			
		END1CYCLE
		
		/*  draw a separating line  */
		if(XA->DrawIt==1)  {
			glColor3fv(CS->Axes);
			glBegin(GL_LINES);
				glVertex2f(0.0f,y-2.5f);
				glVertex2f((GLfloat)gC->Dat->L-(1.0f-xs),y-2.5f);
			glEnd();
		}
			
	END1CYCLE
	
	gsHANDLE_EXTREMA(0.0f,(GLfloat)gC->Dat->L,0.0f, 3.0f * (GLfloat)gC->Dat->M);
	
	
}









void gfnhDrawCategoryProbs(void)
{
	int i,g;
	GLfloat xs = .5f, ys = .5f;
	int n;
	GLfloat w,h;
	GLfloat x, y;
	int Aves = 0;
	GLfloat cum;
	int CW = glutGetWindow();
	char temp[200];
	GLfloat WoH = gfmWoH();

	if(gWindowsSettings[CW]->GFNH_CURRENT_OR_AVES==GFNH_AVES)
		Aves = 1;
	
	gfmSetParsForBarGraphs(gC->Dat->M, gWindowsSettings[CW]->NumCols, xs, ys, &n, &w, &h);
	
	gWindowsSettings[CW]->SelectedElement = gfmBarGraphIdxFromXY(n, gC->Dat->M,  w,  h,  xs,  ys,
					 gWindowsSettings[CW]->XCoord,  gWindowsSettings[CW]->YCoord);
	
	CYCLE_i(gC->Dat)
		cum = 0.0f;
		x = gfmBarGraphXCoord(n,w,xs,i);
		y = gfmBarGraphYCoord(n,h,ys,i);
		
		/*  draw a text label just to the left of the bar, too, with height equal to height of */
		/*  the bar, and moved just left of it (by half of a space character */
		sprintf(temp,"");  /* initialize temp to an empty string */ 
		
		if(gC->Lat->Ind[i]->Flag[Z_FIXED]==1)  {  /* if the individual is from a known category, then make its
												  text label of that same color.  And precede its text label with a "z." */
			glColor3fv(gWindowsSettings[CW]->ColorScheme->Series[gC->Lat->Ind[i]->Z->v /*Flag[FIXED_Z] */]);
			sprintf(temp,"z.");
		}
		else  { /* otherwise make the text label in the standard text color */
			glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		}
		
		/* add an "s." to the text string if it is left out of the mixture sample */
		if(gC->Lat->Ind[i]->Flag[DONT_CONTRIBUTE_TO_PI]==1)  {
			sprintf(temp,"%ss.",temp);
		}
		
		/* concatenate the number on the individual onto its flags */
		sprintf(temp,"%s%d ",temp,i+1);
		
		gfmStrokeString(temp, h ,x, y-h, 3, 0.0f,WoH);
		
		CYCLE_g(gC->Lat)
			
			gfmDrawRectLeaf6f(x,y,w,h, cum, 1.0f, gWindowsSettings[CW]->ColorScheme->Series[g]);
			
			
			if(Aves)
				cum += (GLfloat)gC->Lat->Ind[i]->PofZ[g]->Ave;
			else
				cum += (GLfloat)gC->Lat->Ind[i]->PofZ[g]->v;
		
		END1CYCLE
		
	END1CYCLE
	
	if(gWindowsSettings[CW]->SelectedElement >= 0) {
		sprintf(temp,"Selected Individual = %d", gWindowsSettings[CW]->SelectedElement+1);
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		gfmRenderBitmapString(-.05f,1.1f,GLUT_BITMAP_HELVETICA_12,temp);
	}
}


void gfnhDrawCategoryUniPriProbs(void)
{
	int i,g;
	GLfloat xs = .5f, ys = .5f;
	int n;
	GLfloat w,h;
	GLfloat x, y;
	int Aves = 0;
	GLfloat cum;
	int CW = glutGetWindow();
	char temp[200];
	GLfloat WoH = gfmWoH();

	if(gWindowsSettings[CW]->GFNH_CURRENT_OR_AVES==GFNH_AVES)
		Aves = 1;
	
	gfmSetParsForBarGraphs(gC->Dat->M, gWindowsSettings[CW]->NumCols, xs, ys, &n, &w, &h);
	
	gWindowsSettings[CW]->SelectedElement = gfmBarGraphIdxFromXY(n, gC->Dat->M,  w,  h,  xs,  ys,
					 gWindowsSettings[CW]->XCoord,  gWindowsSettings[CW]->YCoord);
	
	CYCLE_i(gC->Dat)
		cum = 0.0f;
		x = gfmBarGraphXCoord(n,w,xs,i);
		y = gfmBarGraphYCoord(n,h,ys,i);
		
		/*  draw a text label just to the left of the bar, too, with height equal to height of */
		/*  the bar, and moved just left of it (by half of a space character */
		sprintf(temp,"");  /* initialize temp to an empty string */ 
		
		if(gC->Lat->Ind[i]->Flag[Z_FIXED]==1)  {  /* if the individual is from a known category, then make its
												  text label of that same color.  And precede its text label with a "z." */
			glColor3fv(gWindowsSettings[CW]->ColorScheme->Series[gC->Lat->Ind[i]->Z->v /*Flag[FIXED_Z] */]);
			sprintf(temp,"z.");
		}
		else  { /* otherwise make the text label in the standard text color */
			glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		}
		
		/* add an "s." to the text string if it is left out of the mixture sample */
		if(gC->Lat->Ind[i]->Flag[DONT_CONTRIBUTE_TO_PI]==1)  {
			sprintf(temp,"%ss.",temp);
		}
		
		/* concatenate the number on the individual onto its flags */
		sprintf(temp,"%s%d ",temp,i+1);
		
		gfmStrokeString(temp, h ,x, y-h, 3, 0.0f,WoH);
		
		CYCLE_g(gC->Lat)
			
			gfmDrawRectLeaf6f(x,y,w,h, cum, 1.0f, gWindowsSettings[CW]->ColorScheme->Series[g]);
			
			
			if(Aves)
				cum += (GLfloat)gC->Lat->Ind[i]->UniPriPofZ[g]->Ave;
			else
				cum += (GLfloat)gC->Lat->Ind[i]->UniPriPofZ[g]->v;
		
		END1CYCLE
		
	END1CYCLE
	
	if(gWindowsSettings[CW]->SelectedElement >= 0) {
		sprintf(temp,"Selected Individual = %d", gWindowsSettings[CW]->SelectedElement+1);
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		gfmRenderBitmapString(-.05f,1.1f,GLUT_BITMAP_HELVETICA_12,temp);
	}
}





void gfnhDrawPritchardQs(void)
{
	int i,s;
	GLfloat xs = .5f, ys = .5f;
	int n;
	GLfloat w,h;
	GLfloat x, y;
	int Aves = 0;
	GLfloat cum;
	int CW = glutGetWindow();
	char temp[200];
	GLfloat WoH = gfmWoH();

	if(gWindowsSettings[CW]->GFNH_CURRENT_OR_AVES==GFNH_AVES)
		Aves = 1;
	
	gfmSetParsForBarGraphs(gC->Dat->M, gWindowsSettings[CW]->NumCols, xs, ys, &n, &w, &h);
	
	gWindowsSettings[CW]->SelectedElement = gfmBarGraphIdxFromXY(n, gC->Dat->M,  w,  h,  xs,  ys,
					 gWindowsSettings[CW]->XCoord,  gWindowsSettings[CW]->YCoord);
	
	CYCLE_i(gC->Dat)
		cum = 0.0f;
		x = gfmBarGraphXCoord(n,w,xs,i);
		y = gfmBarGraphYCoord(n,h,ys,i);
		
		/*  draw a text label just to the left of the bar, too, with height equal to height of */
		/*  the bar, and moved just left of it (by half of a space character */
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		sprintf(temp,"%d ",i+1);
		gfmStrokeString(temp, h ,x, y-h, 3, 0.0f,WoH);
		
		CYCLE_s(gC->Lat->PritLat)
			
			gfmDrawRectLeaf6f(x,y,w,h, cum, 1.0f, gWindowsSettings[CW]->ColorScheme->Series[s]);
			
			
			if(Aves)
				cum += (GLfloat)gC->Lat->Ind[i]->PritInd->Q[s]->Ave;
			else
				cum += (GLfloat)gC->Lat->Ind[i]->PritInd->Q[s]->v;
		
		END1CYCLE
	END1CYCLE
	
	if(gWindowsSettings[CW]->SelectedElement >= 0) {
		sprintf(temp,"Selected Individual = %d", gWindowsSettings[CW]->SelectedElement+1);
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		gfmRenderBitmapString(-.05f,1.1f,GLUT_BITMAP_HELVETICA_12,temp);
	}
}




void gfnhDrawAlleleFreqs(void)
{
	int k,l;
	GLfloat xs = 1.0f, ys = .5f;
	int n;
	GLfloat w,h,halfh;
	GLfloat x, y;
	int Aves = 0;
	GLfloat cum0,cum1;
	int CW = glutGetWindow();
	char temp[200];
	GLfloat WoH = gfmWoH();
	int MC;

	if(gWindowsSettings[CW]->GFNH_CURRENT_OR_AVES==GFNH_AVES)
		Aves = 1;
		
	/* record the max number of colors in the color scheme */
	MC = gWindowsSettings[CW]->ColorScheme->N;
	
	gfmSetParsForBarGraphs(gC->Dat->L, gWindowsSettings[CW]->NumCols, xs, ys, &n, &w, &h);
	
	halfh = .5f * h;
	
	gWindowsSettings[CW]->SelectedElement = gfmBarGraphIdxFromXY(n, gC->Dat->L,  w,  h,  xs,  ys,
					 gWindowsSettings[CW]->XCoord,  gWindowsSettings[CW]->YCoord);
	
	CYCLE_l(gC->Dat)
		cum0 = 0.0f;
		cum1 = 0.0f;
		
		x = gfmBarGraphXCoord(n,w,xs,l);
		y = gfmBarGraphYCoord(n,h,ys,l);
		
		/*  draw the locus name to the left of the locus */
		sprintf(temp,"%s ",gC->Dat->LocNames[l]);
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		if(l == gWindowsSettings[CW]->SelectedElement)
			glColor3fv(gWindowsSettings[CW]->ColorScheme->Series[0]);  /*  make the selected one a different color */
		gfmStrokeString(temp,h,x,y-h,3,0.0f,WoH);
		
		
		CYCLE_k(gC->Dat)
			
			gfmDrawRectLeaf6f(x,y,w,halfh, cum0, 1.0f, gWindowsSettings[CW]->ColorScheme->Series[k%MC]);
			gfmDrawRectLeaf6f(x,y-halfh,w,halfh, cum1, 1.0f, gWindowsSettings[CW]->ColorScheme->Series[k%MC]);
			
			if(Aves)  {
				cum0 += gC->Lat->Theta[0][l][k]->Ave;
				cum1 += gC->Lat->Theta[1][l][k]->Ave;
			}
			else  {
				cum0 += gC->Lat->Theta[0][l][k]->v;
				cum1 += gC->Lat->Theta[1][l][k]->v;
			}
			
			
		
		END1CYCLE
	END1CYCLE
	
	if(gWindowsSettings[CW]->SelectedElement >= 0) {
		sprintf(temp,"Selected Locus = %s", gC->Dat->LocNames[gWindowsSettings[CW]->SelectedElement]);
		glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
		gfmRenderBitmapString(-.05f,1.1f,GLUT_BITMAP_HELVETICA_12,temp);
	}
}







void gfnhDrawCatProbsHist(void)
{
	int CW = glutGetWindow();
	dval **temp;
	int g;
	char Str[100];
	clipping_struct *Clips = gWindowsSettings[CW]->Clips;
	
	int DE = gWindowsSettings[CW]->DisplayedElement;
	if(DE > gC->Dat->M - 1)  {
		DE = gC->Dat->M - 1;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	if(DE < 0)  {
		DE = 0;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	
	temp = (dval **)ECA_CALLOC((size_t)gC->Dat->Gn->v, sizeof(dval *));
	
	sprintf(Str,"Displayed Individual = %d",DE+1);
	glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
	gfmRenderBitmapString(Clips->xlo + .25f * (Clips->xhi-Clips->xlo),
					Clips->yhi - .22f * (Clips->yhi-Clips->ylo)
					 ,GLUT_BITMAP_HELVETICA_18,Str);
					 
	CYCLE_g(gC->Dat)
		temp[g] = gC->Lat->Ind[DE]->PofZ[g];
	END1CYCLE
	
	gfnhProbsHistogram(&temp,gC->Dat->Gn->v,1);
	
	/*  then draw the axes */
	gfmDrawXYAxes();
	
	free(temp);
}



void gfnhDrawCatUniPriProbsHist(void)
{
	int CW = glutGetWindow();
	dval **temp;
	int g;
	char Str[100];
	clipping_struct *Clips = gWindowsSettings[CW]->Clips;
	
	int DE = gWindowsSettings[CW]->DisplayedElement;
	if(DE > gC->Dat->M - 1)  {
		DE = gC->Dat->M - 1;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	if(DE < 0)  {
		DE = 0;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	
	temp = (dval **)ECA_CALLOC((size_t)gC->Dat->Gn->v, sizeof(dval *));
	
	sprintf(Str,"Displayed Individual = %d",DE+1);
	glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
	gfmRenderBitmapString(Clips->xlo + .25f * (Clips->xhi-Clips->xlo),
					Clips->yhi - .22f * (Clips->yhi-Clips->ylo)
					 ,GLUT_BITMAP_HELVETICA_18,Str);
					 
	CYCLE_g(gC->Dat)
		temp[g] = gC->Lat->Ind[DE]->UniPriPofZ[g];
	END1CYCLE
	
	gfnhProbsHistogram(&temp,gC->Dat->Gn->v,1);
	
	/*  then draw the axes */
	gfmDrawXYAxes();
	
	free(temp);
}


void gfnhDrawAlleleFreqsHist(void)
{
	int CW = glutGetWindow();
	dval ***temp;
	int k,l,a;
	char Str[100];
	clipping_struct *Clips = gWindowsSettings[CW]->Clips;
	
	int DE = gWindowsSettings[CW]->DisplayedElement;
	
	if(DE > gC->Dat->L - 1)  {
		DE = gC->Dat->L - 1;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	if(DE < 0)  {
		DE = 0;
		gWindowsSettings[CW]->DisplayedElement = DE;
	}
	
	/*  set the number of key color items */
	/*  this is really klugie, as it is now.  SHould fix it later. */
	gWinDefs[5]->NumColorKeys = &(gC->Dat->Kl[DE]);
	
	/*  set l, the locus subscript, to the displayed element */
	l=DE;
	
	temp = (dval ***)ECA_CALLOC((size_t)2, sizeof(dval **));
	temp[0] = (dval **)ECA_CALLOC((size_t)gC->Dat->Kl[l], sizeof(dval *));
	temp[1] = (dval **)ECA_CALLOC((size_t)gC->Dat->Kl[l], sizeof(dval *));
	
	sprintf(Str,"Displayed Locus = %d",DE+1);
	glColor3fv(gWindowsSettings[CW]->ColorScheme->Text);
	gfmRenderBitmapString(Clips->xlo + .25f * (Clips->xhi-Clips->xlo),
					Clips->yhi - .22f * (Clips->yhi-Clips->ylo)
					 ,GLUT_BITMAP_HELVETICA_18,Str);
	
	
	CYCLE_k(gC->Dat)
		CYCLE_a
		temp[a][k] = gC->Lat->Theta[a][DE][k];
		END1CYCLE
	END1CYCLE
	
	gfnhProbsHistogram(temp,gC->Dat->Kl[DE],2);
	
	/*  then draw the axes */
	gfmDrawXYAxes();
	
	free(temp[0]);
	free(temp[1]);
	free(temp);
}


/*
	This draws a line curve histogram for probability quantities.
		
	Num1 is the number of histograms in each block of them'
	
	Num2 is the number of blocks of histograms
	
	It uses the hist field in the dval structs, and it also plots the
	current value as a vertical line at the point, under the curve.

*/
void gfnhProbsHistogram(dval ***D, int Num1, int Num2)
{
	int CW = glutGetWindow();
	int i,j;
	GLint factor;
	GLushort pattern;
	
	
	
	/*  cycle through the dvals and plot the hists */
	for(i=0;i<Num2;i++)  { 
		/*  add the stippling here for just i = 0 or i = 1 */
		if(i%2==0) {
			factor = 1;
			pattern = 65535;  /*  no stipple */
		}
		if(i%2==1)  {
			factor = 3;
			pattern = 43690;
		}	
		for(j=0;j<Num1;j++)  {
			gfnhDrawOneHistAsCurveWith2Bars(D[i][j],gWindowsSettings[CW]->ColorScheme->Series[j],factor,pattern,0);
		}
	}
	
}

/*
	draws a line, plus a bar at either end,
	factor and pattern refer to the stippling
	factor = 1 and stipple = 65535 is No stippling.
	AsDensity = 1 means that the heights will be divided by Hist->d so that
	the thing is scaled the way a density would be
*/
void gfnhDrawOneHistAsCurveWith2Bars(dval *D, GLfloat *C, GLint factor, GLushort pattern, int AsDensity)
{
	int i;
	hist_struct *Hist = D->Hist;
	GLfloat d = (GLfloat)D->Hist->d;
	GLfloat x,y,h;
	
	
	/*  set the color */
	glColor3fv(C);
	
	/*  draw rectangles for the bins at the ends (<= Lo and >Hi) */
	/*  draw the one at the left end first */
	if(AsDensity == 1)
		h = (GLfloat)(Hist->H[0]/(Hist->dTot * d));
	else 
		h = (GLfloat)(Hist->H[0]/Hist->dTot);
	glRectf((GLfloat)Hist->Lo, 0.0f, (GLfloat)(Hist->Lo - Hist->d), h );
	
	/*  then the one at the right end */
	if(AsDensity == 1)
		h = (GLfloat)(Hist->H[0]/(Hist->dTot * d));
	else 
		h = (GLfloat)(Hist->H[0]/Hist->dTot);
	glRectf((GLfloat)Hist->Hi, 0.0f, (GLfloat)(Hist->Hi + Hist->d), h );
	
	
	
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(factor,pattern);
	/*  then draw the line */
	glBegin(GL_LINE_STRIP);
		
		for(i=1;i<D->Hist->NumBins-1;i++)  {
			x = d * (GLfloat)( (i-1) + i) / 2.0f;
			y = (GLfloat)(Hist->H[i] / Hist->dTot);
			if(AsDensity==1)
				y /= d;
				
			glVertex2f(x,y);
		}
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	
	/*  then draw the vertical line showing the current value */
	x = (GLfloat)D->v;
	if(D->v >= Hist->Hi)  {
		y = (GLfloat) (Hist->H[Hist->NumBins-1] / Hist->dTot);
	}
	else if( (Hist->HasLeftovers==1) && D->v >= Hist->Left)  {
		y = (GLfloat)(Hist->H[Hist->NumBins-2] / Hist->dTot);
	}
	else if(D->v <= Hist->Lo)  {
		y = (GLfloat)(Hist->H[0] / Hist->dTot);
	}
	else {
		y = (GLfloat)(Hist->H[ (int)((D->v - Hist->Lo) / d) + 1] / Hist->dTot);
	}
	
	if(AsDensity==1)
		y /= d;
		
	glBegin(GL_LINES);
		glVertex2f(x,0);
		glVertex2f(x,y);
	glEnd();

	/*  and that should be it */
}

void gfnhOpenAlleleRightClickWindow(int MenuCall)
{
	int CallingWindow = glutGetWindow();
	int NewWindow;
	
		/*  create the window */
		gfmCreateWindowByIndex(MenuCall);
		
		/*  then update the color scheme and the displayed element in there, */
		/*  as well as whether or not it should draw the legend. */
		NewWindow = glutGetWindow();
		
		gWindowsSettings[NewWindow]->ColorScheme = gWindowsSettings[CallingWindow]->ColorScheme;
		
		gWindowsSettings[NewWindow]->DisplayedElement = gWindowsSettings[CallingWindow]->SelectedElement;
		
		gWindowsSettings[NewWindow]->Legend->DrawLegend = 
				gWindowsSettings[CallingWindow]->Legend->DrawLegend;
	
}


void gfnhOpenCategoryRightClickWindow(int MenuCall)
{
	int CallingWindow = glutGetWindow();
	int NewWindow;
	
	if(MenuCall == 0)  {  /*  this is the call for the histogram window */
		/*  create the window */
		gfmCreateWindowByIndex(6);
		
		/*  then update the color scheme and the displayed element in there, */
		/*  as well as whether or not it should draw the legend. */
		NewWindow = glutGetWindow();
		
		gWindowsSettings[NewWindow]->ColorScheme = gWindowsSettings[CallingWindow]->ColorScheme;
		
		gWindowsSettings[NewWindow]->DisplayedElement = gWindowsSettings[CallingWindow]->SelectedElement;
		
		gWindowsSettings[NewWindow]->Legend->DrawLegend = 
				gWindowsSettings[CallingWindow]->Legend->DrawLegend;
	}
	
}

void gfnhDrawAlleFreqClines(void)
{
	int i,l,ball;
	dval *T;
	char TempStr[1000];
	gsCS;
	gsClips;
	int DE = gWindowsSettings[glutGetWindow()]->DisplayedElement;  /* this will be the locus number */
	double w,h,x,y;  /* width and height of the ball to plot for allele freqs */
	double yvert;
	
	/* first, restrict DE to permissible values */
	if(DE > gC->Dat->L - 1) {
		DE = gC->Dat->L - 1;
	}
	if(DE < 0) {
		DE = -1;
	}
	
	
	/* if DE < 1 then we draw all the averages for all the loci */
	if(DE==-1) {
		for(l=0;l<gC->Dat->L;l++) {  /* cycle over all loci */
			/*  draw the average cline curve */
			glColor3fv(CS->Series[ l % CS->N ]);
			glBegin(GL_LINE_STRIP);
				for(i=0;i<gClines->NumClineXs;i++)  {
					glVertex2f((GLfloat)gClines->ClineX[i],(GLfloat)gClines->Locs[l]->ClineY[i]->Ave);
				}
			glEnd();
		}
		/* then some text */
		glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Text);
		sprintf(TempStr,"Posterior Mean Clines For Each Locus");
		gfmRenderBitmapString(Clips->xlo + .17*(Clips->xhi-Clips->xlo),1.05f,GLUT_BITMAP_HELVETICA_12, TempStr);
	}
	
	else  {  /* otherwise, draw it just for a single locus */
		/* lets draw the "posterior credible set lines"---a bit of a hack, since I just use the observed SD... */
		glColor3fv(CS->Series[2]);
		glBegin(GL_LINE_STRIP);
			for(i=0;i<gClines->NumClineXs;i++)  {
				T = gClines->Locs[DE]->ClineY[i];
				yvert =  QuantileFromHist(T->Hist, T, .05);
				glVertex2f((GLfloat)gClines->ClineX[i],(GLfloat)yvert);
			}
		glEnd();
		glColor3fv(CS->Series[2]);
		glBegin(GL_LINE_STRIP);
			for(i=0;i<gClines->NumClineXs;i++)  {
				T = gClines->Locs[DE]->ClineY[i];
				yvert = QuantileFromHist(T->Hist, T, .95);
				glVertex2f((GLfloat)gClines->ClineX[i],(GLfloat)yvert );
			}
		glEnd();
		

		/* now draw the current cline curve */
			glColor3fv(CS->Series[0]);
			glBegin(GL_LINE_STRIP);
				
				for(i=0;i<gClines->NumClineXs;i++)  {
					glVertex2f((GLfloat)gClines->ClineX[i],(GLfloat)gClines->Locs[DE]->ClineY[i]->v);
		/*			printf("DE = %d   i = %d   x = %f    y = %f\n",DE,i,gClines->ClineX[i],gClines->Locs[DE]->ClineY[i]->v); */
				}
			glEnd();


		/* and on top of that draw the average cline curve */
		glColor3fv(CS->Series[1]);
		glBegin(GL_LINE_STRIP);
			for(i=0;i<gClines->NumClineXs;i++)  {
				glVertex2f((GLfloat)gClines->ClineX[i],(GLfloat)gClines->Locs[DE]->ClineY[i]->Ave);
			}
		glEnd();
		
		/* then draw some text */
		glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Text);
		sprintf(TempStr,"Locus %d.  AveCenter= %.2f  AveBeta= %.5f  Logl = %.4f",DE,gClines->Locs[DE]->alpha->Ave,
						gClines->Locs[DE]->beta->Ave,gClines->Locs[DE]->LogLike->v);
		gfmRenderBitmapString(Clips->xlo + .17*(Clips->xhi-Clips->xlo),1.05f,GLUT_BITMAP_HELVETICA_12, TempStr);

		/* now, draw the allele freq balls */
		w = (Clips->xhi - Clips->xlo) * .02;
		h = (Clips->yhi - Clips->ylo) * .02;
		
		/* then cycle over the locales */
		for(ball=0,i=0;i<=gClines->MaxLocaleIndex;i++)  {
			if(gClines->HasLocale[i]) {
				glColor3fv(CS->Series[ball+4]);
				x = gClines->x[i];
				y = (double)gClines->Locs[DE]->Y[i]->v / (double)gClines->Locs[DE]->N[i]->v;
				glRectf(x-.5*w,y-.5*h,x+.5*w,y+.5*h);
				
				/* then draw one in on the average spot too */
				y = (double)gClines->Locs[DE]->Y[i]->Ave / (double)gClines->Locs[DE]->N[i]->Ave;
				glRectf(x-.5*w,y-.5*h,x+.5*w,y+.5*h);
				
				ball++;
			}
		}	
	}	
	
	
	/*  then draw the axes */
	gfmDrawXYAxes();
	
}


/* the following closes the huge #ifndef COMPILE_NEW_HYB_WITH_NO_GUI block */
#endif 



int main(int argc, char *argv[])
{	
	char File[100];
	int temp,priorchoice;
	double temp_d;
	enum prior_type PiP = JEFFREYS, ThP = JEFFREYS;
	char gnog;
	int DoAsBurnIn = 0, DoAsReps = 0,i,g;
	time_t StartTime, CurrTime, LastTime, RealRepStartTime;
	cli_opts *CL_Opts=NULL;
	int tempSeed1,tempSeed2;
	
	hyb_data *D;
	hyb_prior *P;
	hyb_chain *C;
	
	gPiFixed = 0;  /* default is not fixed */
	
	if(argc>1) {
		CL_Opts = Get_NewHybs_CLI_Opts(argc,argv);
		
		fprintCL_Probs(stdout, CL_Opts);
	}

	else { 
		#ifdef COMPILE_NEW_HYB_WITH_NO_GUI
			fprintf(stderr,"\n\nThis version of NewHybrids was compiled with COMPILE_NEW_HYB_WITH_NO_GUI defined.  Therefore\n");
			fprintf(stderr,"you must use the command-line interface to run it.  Issue the --help option for more information.\n\n");
			exit(1);
		#endif
	
	    /* if no command line options, then do it interactively */
		/* the very first thing that we want to do is get the current working directory
			and store it in a global variable, because when it enters the GLUT 
			interface it is going to forget about that altogether... */
		sprintf(gPWD,"%s/",getenv("PWD"));
		
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

		
		/* Now Process Cline Info */
		ProcessClineOptions(D, gClines);
		
		
		
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

		
		/*  now we can process the individual options a little bit  */
		ProcessIndivOptions(D);
		
		/* put alleles from indivs known to be in Purebred categories into the 
			"PriorizedAllelesPile" */
		PriorizeAllelesFromFixedZ(D);
		
		AddPriorizeAllelesFromFile(D, NULL);
		
			
				
		printf("\n\nGive me two small integers (>0) for random number seeds\n\n");
		scanf("%d%d",&tempSeed1,&tempSeed2);
		setall((long)tempSeed1,(long)tempSeed2);
		
		
		
		printf("\n\nEnter the choices for prior type for pi (mixing proportion):\n\t0 for Jeffreys\n\t1 for Uniform\n\t2 to specify it with fixed values\n->");
		scanf("%d",&priorchoice);


		switch(priorchoice) {
			case(0): 
				PiP = JEFFREYS;
				break;
			case(1):
				PiP = UNIFORM;
				break;
			case(2):
				PiP=FIXED_PRIOR;  /* I have to make a new enum for this that reports it properly */
				gPiFixed = 1; 
				gPiFixedValues = (double *)ECA_CALLOC(D->Gn->v,sizeof(double));
				temp_d = 0.0;
				printf("\nEnter values for fixed components of Pi.\nThese will be rescaled so as to sum to 1.0 if necessary.\n");
				CYCLE_g(D)
					printf("%s?  --> ",D->CategoryNames[g]);
					scanf("%lf",&(gPiFixedValues[g]));
					temp_d += gPiFixedValues[g];
				END1CYCLE
				
				printf("\nThank You! Using values:\n");
				CYCLE_g(D)
					gPiFixedValues[g] /= temp_d;
					printf("%s:  %f\n",D->CategoryNames[g],gPiFixedValues[g]);
				END1CYCLE
				break;
			default:
				PiP = JEFFREYS;
				printf("\n\nInvalid choice.  Pi Prior set to Jeffreys by default.\n\n");
		}


		P = CreatePriors(D,PiP,ThP);
		C = CreateLatentChain(D,P);
		InitializeChain(C);
		C->Seed1 = tempSeed1;
		C->Seed2 = tempSeed2;
		
		printf("\n\nData all read and ready!!!\n\n");
		
		printf("\n\nEntering GLUT interface...\n\n");
	}
	
	
	if(CL_Opts != NULL)  {  /* if we didn't do it all interactively, then we have to copy some variables, etc., over */
		/*  get the DATA: */
		printf("DATA_INITIALIZATION: Preparing to read data from %s\n",CL_Opts->DataFilePath);
		D = GetData(CL_Opts->DataFilePath);
		
		/* Now Process Cline Info */
		printf("DATA_INITIALIZATION: Processing cline options (if any)\n");
		ProcessClineOptions(D, gClines);
		
		/* then initialize the gtyp frequency category probabilities */
		printf("DATA_INITIALIZATION: Initializing genotype frequency category probabilities\n");
		if(strlen(CL_Opts->GtypCatFilePath)>0) {
			GetGtypFreqCats(D,CL_Opts->GtypCatFilePath);
		}
		else {
			CopyGtypFreqCatsFromCL(D,CL_Opts);
		}
		
		/*  now we can process the individual options a little bit  */
		ProcessIndivOptions(D);
		
		/* put alleles from indivs known to be in Purebred categories into the 
			"PriorizedAllelesPile" */
		PriorizeAllelesFromFixedZ(D);
		
		/* then add them in from a file, if there is one */
		if(strcmp(CL_Opts->AlleFreqPriorPath,"UNSET") != 0) {
			AddPriorizeAllelesFromFile(D, CL_Opts->AlleFreqPriorPath);
		}
		
		/* and copy across information about the pi priors */
		PiP = CL_Opts->PiPriType;
		ThP = CL_Opts->ThetaPriType;
		if(PiP==FIXED_PRIOR) { int i;
			gPiFixed = 1; 
			gPiFixedValues = (double *)ECA_CALLOC(D->Gn->v,sizeof(double));
			for(i=0;i<D->Gn->v;i++)  {
				gPiFixedValues[i] = CL_Opts->PiFixedValues[i];
			}
		}
		
		
		/* copy across the Trace Report request to the D structure */
		D->PiTraceReport = CL_Opts->PiTraceReport;
        D->ZTraceReport = CL_Opts->ZTraceReport;
		
		
		P = CreatePriors(D,PiP,ThP);
		C = CreateLatentChain(D,P);
		
		if(CL_Opts->Seed1==0 && CL_Opts->Seed2==0) {
			SeedFromFile("newhyb_seeds");
		}
		else {
			printf("SEEDS_FROM_COMMAND_LINE : %ld  %ld\n",CL_Opts->Seed1,CL_Opts->Seed2);
			setall(CL_Opts->Seed1,CL_Opts->Seed1);
			C->Seed1 = (int)CL_Opts->Seed1;
			C->Seed2 = (int)CL_Opts->Seed2;
		}
		
		InitializeChain(C);
	}
	

	/* Here this is set up so that if you use the interactive version, you end up having to use the
	GUI, for now.  You have the choice if you use the command line */
	if(CL_Opts==NULL || CL_Opts->NoGui==0) {
	
		#ifndef COMPILE_NEW_HYB_WITH_NO_GUI
		/*  make the global pointer point to where it ought to */
		gC = C;
		/*  then allocate space for a sliding trace struct that holds up to 1,000 entries */
		gCompleteDataLogLike = gfduAllocSlidingTrace(1000,1);
			

		/*  here must call glutInit from within main, with the command line option */
		/*  pointer, and then call gfmInitGFM, and that does it.   */
		glutInit(&argc, argv);
		gfmInitGFM();
		#else
			fprintf(stderr,"\n\nThis version of NewHybrids was compiled with COMPILE_NEW_HYB_WITH_NO_GUI defined.  Therefore\n");
			fprintf(stderr,"it MUST be invoked with the --no-gui option.  Try adding --no-gui to your last command line.  Cheers.\n\n");
			exit(1);
		#endif
	}
	else {int dummy;
		dummy = RunWithoutGraphics(C, CL_Opts->NumBurnIn, CL_Opts->NumPostBurnIn);
	}


	return(1);
}
