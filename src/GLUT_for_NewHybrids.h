/*  some #defines to make it simple to use the integer options */
#define GFNH_CURRENT_OR_AVES IntegerOptions[0]
#define GFNH_CURRENT 0
#define GFNH_AVES	 1

#define GFNH_MAX_CATEGORIES 100


/*  global pointer to the data structure */
hyb_chain *gC;
/* a global variable to hold the sliding trace struct for the complete data loglikelihood */
GLOB sliding_trace_struct *gCompleteDataLogLike;

/*  function prototypes */
void gfnhDrawConsoleWindow(void);
void gfnhConsoleKeys(unsigned char key, int x, int y);
void gfnhDrawCategoryProbs(void);
void gfnhDrawAlleleFreqs(void);
void gfnhProbsHistogram(dval ***H, int Num1, int Num2);
void gfnhDrawCatProbsHist(void);
void gfnhDrawOneHistAsCurveWith2Bars(dval *D, GLfloat *C, GLint factor, GLushort pattern, int AsDensity);
void gfnhDrawAlleleFreqsHist(void);
void gfnhDrawObservedData(void);
void gfnhDrawWs(void);
void gfnhDrawW_Aves(void);
void gfnhDrawCompleteDataLogLTrace(void);

void gfnhDrawCategoryUniPriProbs(void);
void gfnhDrawCatUniPriProbsHist(void);

void gfnhOpenAlleleRightClickWindow(int MenuCall);
void gfnhOpenCategoryRightClickWindow(int MenuCall);
void gfnhDrawKullbLeib(void);

void gfnhDrawPritchardQs(void);
void gfnhDrawPiHists(void);
void gfnhDrawPiValues(void);



void gfnhDrawAlleFreqClines(void);