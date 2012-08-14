/*  some #defines to make it simple to use the integer options */
#define GFNH_CURRENT_OR_AVES IntegerOptions[0]
#define GFNH_CURRENT 0
#define GFNH_AVES	 1

#define GFNH_MAX_CATEGORIES 100


/*  for controlling the type of simulation we are doing. */
/*  this will control things in PritchSingleSweep.  It will  */
/*  not have any effect on variable initialization and memory allocation... */
/*  That will all still be done as if we are doing it JABES. */
/*  choices for this are defined in NewHybrids.h: */
/*  PRITCHARD and JABES */
GLOB int TypeOfPritchEtcSim;




/*  global pointer to the data structure */
GLOB hyb_chain *gC;


/*  function prototypes */
void gfnhDrawConsoleWindow(void);
void gfnhDrawCategoryProbs(void);
void gfnhDrawAlleleFreqs(void);
void gfnhProbsHistogram(dval ***H, int Num1, int Num2);
void gfnhDrawOneHistAsCurveWith2Bars(dval *D, GLfloat *C, GLint factor, GLushort pattern, int AsDensity);
void gfnhDrawAlleleFreqsHist(void);
void gfnhDrawObservedData(void);
void gfnhDrawCompleteDataLogLTrace(void);

void gfnhOpenAlleleRightClickWindow(int MenuCall);
void gfnhOpenPritchQRightClickWindow(int MenuCall);
void gfnhDrawKullbLeib(void);

void gfnhDrawPritchardQs(void);



void gfpeDrawAlphaHists(void);
void gfpeDrawPritchQHist(void);
void gfpeDrawXiHists(void);
void gfpeDrawPiHists(void);

void gfpehDrawPofV(void);