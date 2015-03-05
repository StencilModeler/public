#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "cycle.h"
#include "platforms.h"
#include "common.h"
#include "opts.h"

#ifdef HAVE_PAPI
#include <papi.h>
#endif

/*
 * TODO: explain parameters
 * Minimum time that could spent to bring data to CPU: Optimistic
 *
 * TcoldRAM: Compute 1 plane (size II*JJ) using non-prefetched data
 *           This means ReadPlanes planes using TclineRAMrNP
 * TstreamRAM: Compute K-1 remaining planes (II*JJ) using prefetched data
 *             This means only 1 plane for each new plane to compute
 *             using TclineRAMr because remaining data is already in cache.
 * TwriteRAM: If the platform uses Write-back (write-allocate) policy,
 *            data must be brought from RAM to cache, before saving it.
 */
double transfer_optimistic( int nx, int ny, int nz, int length, int timesteps )
{
  /* General */
  int I, J, K;
  int II, JJ;

  /* CPU - RAM */
  double TcoldRAM, TstreamRAM, TwriteRAM;
  double TclineRAMr, TclineRAMrNP, TclineRAMw;
  double TtotalRAM;
  double Ttotal;

  DEBUG_FPRINTF( stderr, "OPTIMISTIC (CASE 1): ######################################\n" );

  I = nx;
  J = ny;
  K = nz;
  II = I+2*length;
  JJ = J+2*length;

  /* RAM */
  TclineRAMr   = cacheline / bwRAMr;
  TclineRAMrNP = cacheline / bwRAMrNP;
  TclineRAMw   = cacheline / bwRAMw;
  TcoldRAM     = (ReadPlanes * ceil(II / elemCLine) * JJ) * TclineRAMrNP;
  TstreamRAM   = ceil(II / elemCLine) * JJ * (K-1) * TclineRAMr;
#if WRITEBACK
  TwriteRAM  = WritePlanes * ceil(I / elemCLine) * J * K * TclineRAMw;
#else
  TwriteRAM  = 0.0;
#endif
  TtotalRAM  = TcoldRAM + TstreamRAM + TwriteRAM;

  Ttotal = TtotalRAM;
  DEBUG_FPRINTF( stderr, "Ttotal \tTcoldRAM \tTstreamRAM \tTwriteRAM\n");
  DEBUG_FPRINTF( stderr, "%g \t%g \t%g \t%g\n", Ttotal, TcoldRAM, TstreamRAM, TwriteRAM);

  DEBUG_FPRINTF( stderr, "OPTIMISTIC (CASE 1): ######################################\n\n" );


  return Ttotal;
}

/*
 * TODO: explain parameters
 * Maximum time that could spent to bring data to CPU: Pesimistic
 *
 * Tcold: Bring data from RAM using CASE4 (2*nplanes-1) planes
 *        to compute an output plane (non-prefetch bandwidth).
 * TwriteRAM: If the platform uses Write-back (write-allocate) policy,
 *            data must be brought from RAM to cache, before saving it.
 */
double transfer_pesimistic( int nx, int ny, int nz, int length, int timesteps )
{
  /* General */
  int I, J, K;
  int II, JJ, KK;

  /* CPU - RAM */
  double TcoldRAM, TstreamRAM, TwriteRAM;
  double TclineRAMrNP, TclineRAMw;
  double TtotalRAM;
  double Ttotal;

  DEBUG_FPRINTF( stderr, "PESIMISTIC (CASE 2): ######################################\n" );

  I = nx;
  J = ny;
  K = nz;
  II = I+2*length;
  JJ = J+2*length;
  KK = K+2*length;

  /* RAM */
  TclineRAMrNP = cacheline / bwRAMrNP;
  TclineRAMw   = cacheline / bwRAMw;
  TcoldRAM     = ((2 * ReadPlanes - 1) * ceil(II / elemCLine) * JJ) * KK * TclineRAMrNP;
  TstreamRAM   = 0.0F;
#if WRITEBACK
  TwriteRAM  = WritePlanes * ceil(I / elemCLine) * J * K * TclineRAMw;
#else
  TwriteRAM  = 0.0F;
#endif
  TtotalRAM  = TcoldRAM + TstreamRAM + TwriteRAM;

  Ttotal = TtotalRAM;
  DEBUG_FPRINTF( stderr, "Ttotal \tTcoldRAM \tTstreamRAM \tTwriteRAM\n");
  DEBUG_FPRINTF( stderr, "%g \t%g \t%g \t%g\n", Ttotal, TcoldRAM, TstreamRAM, TwriteRAM);

  DEBUG_FPRINTF( stderr, "PESIMISTIC (CASE 2): ######################################\n\n" );

  return Ttotal;
}

/*
 * Interpolates a value (y) given a range of data (lowerBound, upperBound) for
 * x and y axes:
 *
 *  x -> E [x0, x1]
 *  y -> unknown E [y0, y1]
 *
 * Methods:
 *  * No interpolation (truncation)
 *  * Linear interpolation
 *  * Exponential interpolation
 *  * Logarithmic interpolation
 *
 * Return y as the interpolated data [y0, y1]
 */
double interpolate( int method, double x, double x0, double x1, double y0, double y1 )
{
  /* Interpolated value */
  double y;
  double k;


  /* Interpolate given a method (no interp, linear, etc...) */
  switch( method ) {

    /* No interpolation at all - return lowerBound */
    case NO_INTERP:
      y = y0;
      break;

    /* Linear interpolation */
    case LIN_INTERP:
      y = y0 + ((y1 - y0) * ((x - x0) / (x1 - x0)));
      break;

    /* Exponential interpolation */
    /* y = (e^kx - 1) / (e^k - 1)   k E (0,+inf)
     * or
     * y = x / (1 + a * (1 - x))    a E [0,+inf) */
    case EXP_INTERP:
#ifndef KINTERP
      k = exp(1);
#else
      k = KINTERP;
#endif
      y = y0 + ((y1 - y0) * ((exp(k*((x-x0)/(x1-x0))) - 1) / (exp(k) - 1)));
      break;

    /* Logarithmic interpolation */
    /* y = ( 1 - e^(-kx) ) / ( 1 - e^(-k) )   k E (0,+inf) */
    case LOG_INTERP:
#ifndef KINTERP
      k = exp(1);
      //k = 1; //For case 4
#else
      k = KINTERP;
#endif
      y = y0 + ((y1 - y0) * ((1 - exp(-k*((x-x0)/(x1-x0)))) / (1 - exp(-k))));
      break;

    default:
      fprintf( stderr, "Error: interpolation method not supported!\n" );
      exit(EXIT_FAILURE);
  }

  DEBUG_FPRINTF(stderr, "***** Interpolation: (method: %d)\n", method);
  DEBUG_FPRINTF(stderr, "\ty: %g, for x: %g, x0: %g, x1: %g, y0: %g, y1: %g\n",
                y, x, x0, x1, y0, y1);

  /* Return interpolated data */
  return y;
}

/* Returns the number of misses by means of number of planes to read
 * for each plane processed of the stencil computation */
double read_misses( char *textLx, int I, int II, int J, int K, int length,
                    double sizeLx, double assocLx, int interp, int *Kread, double *writePlanes, double *nplanesLxNCI )
{
  int JJ = J+2*length;
  int KK = K+2*length;

  double nplanes, nplanesSize;
  double planeSizeW = I * J;   // Plane size used for writing
  double planeSizeR = II * JJ; // Plane size used for reading

  double assocSize;
  double nplanesLx;

  int lowerRule, upperRule;
  double x, x0, x1, y0, y1;
  double scaleII, scaleJJ;

  /* Rules */
  int R1, R2, R3, R4, R5;


  /* TODO: Don't consider planeSizeW if write-through policy */
#if WRITEBACK
  nplanes     = ReadPlanes + WritePlanes;                           // Total number of planes keept in caches
  nplanesSize = ReadPlanes * planeSizeR + WritePlanes * planeSizeW; // Total kB of all planes in memory

  /* By default writePlanes is the macro */
  *writePlanes = WritePlanes;
#else
  nplanes     = ReadPlanes;
  nplanesSize = ReadPlanes * planeSizeR;
#endif
  DEBUG_FPRINTF( stderr,"Total memory (elems) per stencil plane: %g (%.0f planes)\n",
                 nplanesSize, nplanes );


  /* Associativies that every plane requires for non fitting cases */
  /* How much space/assoc requires the central plane */
  /* assoc = MIN(planeSizeR/((sizeL1/word)/assocL1), assocL1); */   // Number of associativities required
  assocSize  = (sizeLx/word)/assocLx;                 // Size per associativity
  DEBUG_FPRINTF( stderr,"\t%s assocSize: %g\n", textLx, assocSize );


/* Rules definition */
#if INTERP == NO_INTERP
  R1 = ((sizeLx/word) >= nplanesSize);
  R2 = 0; //FALSE only for interpolation case;
  R3 = (((sizeLx/word)*COLUMN_RATIO) >= planeSizeR);
  R4 = (((sizeLx/word)-assocSize) >= ReadPlanes * II);
#else
  //if (((sizeLx/word) - assocSize) >= nplanesSize) { //if ((sizeLx/word) >= nplanesSize) {
  R1 = ((sizeLx/word)*COLUMN_RATIO >= nplanesSize);
  R2 = ((sizeLx/word) >= nplanesSize);
  R3 = (((sizeLx/word)*COLUMN_RATIO) >= planeSizeR);
  //if (((sizeLx/word)-assocSize) < ReadPlanes * II) { //if (((sizeLx/word)) < ReadPlanes * II) {
  R4 = (((sizeLx/word)*COLUMN_RATIO) >= ReadPlanes * II);
  //R4 = ((((sizeLx/word)-assocSize)*COLUMN_RATIO) >= ReadPlanes * II);
#endif

  /* By default K planes to read */
  *Kread = K;

  /* If all the planes required fit in cache (ReadPlanes + WritePlanes (1)):                        */
  /* Total reuse by planes: Fetch 1 plane for reading + 1 plane for writing (write allocate policy) */
  if (R1) {
    /* At least only nplanesLx = 1 of misses, but coef, indexes could push up this value a little more for L1 */
    nplanesLx = 1.0;
    DEBUG_FPRINTF( stderr,"\t%s: Total reuse (C1), nplanes: %f\n", textLx, nplanesLx);
    *Kread = KK;

    lowerRule = 0;
    upperRule = 1;
    goto endRules;
  }

  /* Transition between 1 (R1) -> (Pread - 1) (R3): Interpolation case */
  if (R2) {
    /* Initial number of planes to read if total data does not fit in cache */
    nplanesLx = ReadPlanes - 1.0;
    DEBUG_FPRINTF( stderr,"\t%s: No total reuse on K, reuse on central plane (C1 U C2), nplanes: %f\n", textLx, nplanesLx);

    lowerRule = 1;
    upperRule = 2;
    goto endRules;
  }

#if SEMI
  /* Semi-stencil: writePlanes is increased to 3 if central plane does not fit in cache */
  *writePlanes+=1.0;
  DEBUG_FPRINTF( stderr,"\t%s: writePlanes: %f\n", textLx, *writePlanes);
#endif

  /* If central plane fits in half of the cache (COLUMN_RATIO: columns central/columns in total) */
  /* the number of planes should be (ReadPlanes - 1) because reuse central plane                 */
  if (R3) {
    /* Miss back and front planes only */
    nplanesLx = ReadPlanes - 1.0;
    DEBUG_FPRINTF( stderr,"\t%s: No reuse on K, reuse on central plane (C2 U C3), nplanes: %f\n", textLx, nplanesLx);

    lowerRule = 2;
    upperRule = 3;
    goto endRules;
  }

  /* If central plane does not fit in half of the cache (COLUMN_RATIO: columns central/columns in total) */
  /* the number of planes should be ReadPlanes                                                           */
  if (R4) {
    nplanesLx = ReadPlanes;
    DEBUG_FPRINTF( stderr,"\t%s: No reuse on K (neither central plane) (C3 U C4), nplanes: %f\n", textLx, nplanesLx);

    lowerRule = 3;
    upperRule = 4;
    goto endRules;
  }

  /* If ncolumns of central plane do not fit in cache * FACTOR */
  /* ReadPlanes extras has to be read for each iteration on k  */
  /* They should be reused in J iteration                      */
  /* Miss all planes (included central plane: no reuse columns on central plane) */
  nplanesLx = (2.0 * ReadPlanes) - 1.0;
  DEBUG_FPRINTF( stderr,"\t%s: No reuse, miss all (planes, central plane and columns) (C4), nplanes: %f\n", textLx, nplanesLx);

  lowerRule = 4;
  upperRule = 5;
endRules:

  /* Interpolation cases */
  if ( interp != NO_INTERP ) {
    /* Compute upper and lower bound for interpolation - Depends on rules */
    switch( lowerRule ) {

    /* nplanesLx = 1 */
    case 0:

      /* No interpolation for Rule 1 */
      x0 = x = 1.0;
      x1 = 2.0;

      /* Lower and upper bounds for Y */
      y0 = 1.0;  /* 1 */
      y1 = 1.0;  /* 1 */
      break;

    /* nplanesLx = 1 -> upwards (P-1): C1 U C2 */
    case 1:

      /* ((sizeL1/word)*COLUMN_RATIO <= nplanesSize (Rule 1) */
      //x0 = (((sizeLx/word) - assocSize) - (I * J)) / (ReadPlanes * JJ);
      x0 = (((sizeLx/word)*COLUMN_RATIO) - ((*writePlanes) * I * J)) / (ReadPlanes * JJ);

#if 1
      /* ((sizeLx/word) * COLUMN_RATIO) < planeSizeR   (Rule 3) */
      x1 = (((sizeLx/word)*COLUMN_RATIO) / JJ); // Old rule
#else
      /* ((sizeLx/word) >= nplanesSize) (Rule 2) */
      x1 = ((sizeLx/word) - ((*writePlanes) * I * J)) / (ReadPlanes * JJ);
#endif

      /* Current point to interpolate */
      x = MAX(x0,MIN(II,x1));

      /* Lower and upper bounds for Y */
      y0 = 1.0;               /* 1 */
      y1 = ReadPlanes - 1.0;  /* P - 1 */
      break;

    /* nplanesLx = P-1 -> upwards (P): C2 U C3 */
    case 2:

      /* ((sizeL1/word) <= nplanesSize (Rule 2) */
      //x0 = (((sizeLx/word) - assocSize) - (I * J)) / (ReadPlanes * JJ);
      x0 = ((sizeLx/word) - ((*writePlanes) * I * J)) / (ReadPlanes * JJ);

      /* ((sizeLx/word) * COLUMN_RATIO) < planeSizeR   (Rule 3) */
      x1 = (((sizeLx/word) * COLUMN_RATIO) / JJ);

      /* Current point to interpolate */
      x = MAX(x0,MIN(II,x1));

      /* Lower and upper bounds for Y */
      y0 = ReadPlanes - 1.0;  /* P - 1 */
      y1 = ReadPlanes;        /* P */
      break;

    /* nplanesLx = P: C3 */
    case 3:

      /* No interpolation for Rule 3 */
      x0 = x = 1.0;
      x1 = 2.0;

      /* Lower and upper bounds for Y */
      y0 = ReadPlanes;                /* P */
      y1 = ReadPlanes;                /* P */
      break;

    /* nplanes Lx = P -> upwards (2P - 1): C3 U C4 */
    case 4:

      /* ((sizeLx/word) * COLUMN_RATIO) < planeSizeR (Rule 3) */
      x0 = (((sizeLx/word) * COLUMN_RATIO) / JJ);

      /* ((sizeL1/word) * FACTOR) < ReadPlanes * II (Rule 4) */
      //x1 = (((sizeLx/word) - assocSize) / ReadPlanes);
      x1 = ((sizeLx/word) / ReadPlanes);
      //x1 = (((sizeLx/word)*COLUMN_RATIO) / ReadPlanes);

      /* Current point to interpolate */
      x = MAX(x0,MIN(II,x1));

      /* Lower and upper bounds for Y */
      y0 = ReadPlanes;                /* P */
      y1 = (2.0 * ReadPlanes) - 1.0;  /* 2P - 1 */
      break;

    /* Just in case */
    default:
      y0 = nplanesLx;
    }


    DEBUG_FPRINTF(stderr, "lowerRule = %d, upperRule = %d\n", lowerRule, upperRule);

    DEBUG_FPRINTF(stderr, "%s: (%d,%d,%d), %d interpolation, nplanesLx: %f\n", textLx,
                  I, J, K,  interp, nplanesLx);

    /* Interpolate between rules */
    nplanesLx = interpolate( interp, x, x0, x1, y0, y1 );
  }


  /* There is a factor of conflict and capacity misses given by:
   *
   * Interferences between central planes and remaining depending
   * on JJ factor (J + 2*LENGTH) and II factor
   *
   *    scaleII = (1 - (Pread - 1)/(II)) E (0,1]
   *
   *    scaleJJ = (1 - (Pread - 1)/(JJ)) E (0,1]
   *
   * */
  *nplanesLxNCI = nplanesLx;
#ifndef NO_SCALE
  //scaleII = (1.0 - ((ReadPlanes - 3.0)/II));
  scaleII = 1.0;
  scaleJJ = (1.0 - ((ReadPlanes - 1.0)/JJ));
  //scaleJJ = (1.0 - (1.0/JJ));
  //scaleJJ = (1.0 - ((ReadPlanes)/JJ));
  //scaleJJ = 1.0;
  //scaleJJ = ((II*(JJ-ReadPlanes))/(II*JJ));
  if ((lowerRule > 1) && (upperRule > 1)) nplanesLx *= (scaleII*scaleJJ);
  //else nplanesLx += (scaleII*scaleJJ);
  //if ((lowerRule > 1) && (upperRule > 1)) nplanesLx *= (scaleJJ);
  //else nplanesLx += (scaleII * scaleJJ);
  DEBUG_FPRINTF(stderr, "scaleII = %f, scaleJJ = %f, nplanesLx = %f\n", scaleII, scaleJJ, nplanesLx);
#endif


  /* Return number of planes misses for such configuration */
  return nplanesLx;
}

/*
 * TODO: explain parameters
 * This function computes the estimated time to compute
 * the stencil computation.
 */
double transfer_ideal( int nx, int ny, int nz,
                       int tx, int ty, int tz,
                       int length, int timesteps )
{
  /* General */
  int I, J, K, Kread;
  int II, JJ;
  double NB, NBI, NBJ, NBK;

  double y0, y1;
  double volumeRead, volumeReadCL, volumeWriteCL;

  /* CPU - L1 */
  double regLoopR, regLoopW;
  double regConst, regIndex;
  double TstreamL1, TnstreamL1;
  double TwordL1r, TwordL1w;
  double TwordL1rNP;//, TwordL1wNP;
  double TtotalL1;
  double MissL1P, MissL1NP, HitL1P, HitL1NP; 
  double nplanesL1P, nplanesL1NP;
  double nplanesL1, nplanesL1NCI;

  /* L1 - L2 */
  double TstreamL2, TnstreamL2;
  double TclineL2r;//, TclineL2w;
  double TclineL2rNP; //, TclineL2wNP;
  double TtotalL2;
  double MissL2P, MissL2NP, HitL2P, HitL2NP; 
  double nplanesL2P, nplanesL2NP;
  double nplanesL2, nplanesL2NCI;

  /* L2 - L3 */
  double TstreamL3, TnstreamL3;
  double TclineL3r; //, TclineL3w;
  double TclineL3rNP; //, TclineL3wNP;
  double TtotalL3;
  double MissL3P, MissL3NP, HitL3P, HitL3NP;
  double nplanesL3P, nplanesL3NP;
  double nplanesL3, nplanesL3NCI;

  /* L3 - RAM */
  double TstreamRAM, TnstreamRAM;
  double TclineRAMr, TclineRAMw;
  double TclineRAMrNP, TclineRAMwNP;
  double TtotalRAM;
  double Ttotal;

  /* WRITE */
  double nplanesWP, nplanesWNP;
  double writePlanes = WritePlanes; // For Semi-stencil case on Writing (2->3)
  double MissWP, MissWNP;
  double Twrite, TwriteP, TwriteNP;

  /* Prefetching Efficiency coefficients. Range between [0,1].
   *  - 0 means no prefetching at all.
   *  - 1 means prefetching everything. */
  double pEL1, pEL2, pEL3 = 0.0;
  double prefBlocking;
  double prefEL1[] = prefEffL1;
  double prefEL2[] = prefEffL2;

  DEBUG_FPRINTF( stderr, "IDEAL (CASE 3): ###########################################\n" );

#if 1
/* Compute blocking parameters */
  NBI = nx/(double)tx;
  NBJ = ny/(double)ty;
  NBK = nz/(double)tz;
  NB  = (NBI * NBJ * NBK);
  I = tx;
  J = ty;
  K = tz;


  // If NBI > 1 or VECTOR processor -> words are also fetched multiple of cline
  if ((NBI > 1) && VECTOR) {
    I  = ceil((I+UNALIGNED)/elemCLine)*elemCLine;
    II = ceil((I+2*length)/elemCLine)*elemCLine;
  }
  else
    II = I+2*length;


#if 0
#define TP 2
#define LAP 20
  /* Blocking effect in Prefetched planes */
  if (NB > 1) {
    if (I > TP*elemCLine) { // Trigger Prefetching activation
      I  += LAP*elemCLine; // Look-ahead Prefetching size
      II += LAP*elemCLine; // Look-ahead Prefetching size
      prefBlocking = ((TP+LAP)*elemCLine)/II;
    }
    else prefBlocking = 0.0;  // Prefetching is not triggered
  }
  else prefBlocking = 1.0;
  DEBUG_FPRINTF( stderr, "prefBlocking: %f\n", prefBlocking);
#endif


#if 0
#if (VECTOR == 1)
  if (length < VECSIZE) II = I+2*VECSIZE;
  else                  II = I+2*length;
#endif
#endif

  JJ = J+2*length;

  DEBUG_FPRINTF( stderr, "Blocking is ON\n" );
  DEBUG_FPRINTF( stderr, "II: %d, JJ: %d, NBI: %f, NBJ: %f, NBK: %f, NB: %f\n",
                 II, JJ, NBI, NBJ, NBK, NB);

#else
  /* No blocking */
  NBI = 1;
  NBJ = 1;
  NBK = 1;
  NB  = 1;
  I = nx;
  J = ny;
  K = nz;
  II = I+2*length;
  JJ = J+2*length;
#endif

  /* Non-streamed for blocking */
// if (I != nx) {
//   Tfirst = (nx * ny * nz)/I;
// }
// if ((I == nx) && (ny != J)) {
//   Tfirst = (nx * ny * nz)/(I*J);
// }
// if ((I == nx) && (J == ny)) {
//   Tfirst = (nx * ny)/(I*J);
// }


  /******************************************/
  /* CPU - L1: Cost transfer from L1 to CPU */
  /******************************************/
  DEBUG_FPRINTF( stderr, "*** CPU - L1 ***\n" );
  TwordL1r  = word / bwL1r;
  TwordL1w  = word / bwL1w;
  TwordL1rNP = word / bwL1rNP;
  //TwordL1wNP  = word / bwL1wNP;

  /* Approximate numbers of registers used to compute the internal Loop */
  regLoopR = ReadWords;
  regLoopW = WriteWords;
  regConst = length * dim + 1;          // FPR registers used for constants
  regIndex = 2 * length * (dim-1) + 1;  // GPR registers used for indexes/offsets
  if (regFPR < regConst) {
    regLoopR += regConst - regFPR; // Cannot keep constants in FPR registers
  }
  if (regGPR < regIndex) {
    regLoopR += regIndex - regGPR; // Cannot keep offsets/indexes in GPR registers
    regLoopW += regIndex - regGPR;
  }


  /* Compute Read Misses - L1 */
  nplanesL1 = read_misses( "L1", I, II, J, K, length, sizeL1, assocL1, INTERP, &Kread, &writePlanes, &nplanesL1NCI);

  volumeRead    = II * JJ * Kread * NB;                             // Read volume in points
  volumeReadCL  = ceil((II * JJ * Kread) / elemCLine) * NB;         // Read volume in Cache Lines
#if WRITEBACK
  volumeWriteCL = ceil((I * J * K) / elemCLine) * NB;               // Write volume in Cache Lines (write-allocate policy)
#else
  volumeWriteCL = 0.0;                                              // Write volume (write-through policy)
#endif


  /* Compute prefetching efficiency for L1 */
  /* For architectures with Prefetching system - Store reads may be prefetched as well - WRITEBACK */
  if ((prefL1) && (WRITEBACK)) {
    /* Add writePlanes because Writeback prefetch */
    y0 = prefEL1[MIN(((int)nplanesL1)+((int)writePlanes),   nPrefEffL1)-1];
    y1 = prefEL1[MIN(((int)nplanesL1)+((int)writePlanes)+1, nPrefEffL1)-1];
  }
  else {
    y0 = prefEL1[MIN(((int)nplanesL1),   nPrefEffL1)-1];
    y1 = prefEL1[MIN(((int)nplanesL1)+1, nPrefEffL1)-1];
  }
  pEL1 = interpolate( INTERP, nplanesL1, (int)nplanesL1, ((int)nplanesL1)+1, y0, y1 );


#if 1
  nplanesL1P  = nplanesL1 * pEL1;         // Planes to read WITH prefetching
  nplanesL1NP = nplanesL1 * (1.0 - pEL1); // Planes to read WITHOUT prefetching
#else
  nplanesL1NP = MAX(nplanesL1 - prefL1, 0.0); // Planes to read WITHOUT prefetching
  nplanesL1P  = nplanesL1 - nplanesL1NP;      // Planes to read WITH prefetching
#endif


  MissL1P  = nplanesL1P * volumeReadCL;
  MissL1NP = nplanesL1NP * volumeReadCL;
  if (pEL1 > 0.0) {
    HitL1P   = MAX((regLoopR * I * J * K * NB) - (nplanesL1 * volumeRead), 0.0);
    HitL1NP  = 0;
  }
  else {
    HitL1P   = 0;
    HitL1NP  = MAX((regLoopR * I * J * K * NB) - (nplanesL1 * volumeRead), 0.0);
    //HitL1NP  = (regLoopR * I * J * K * NB);
  }

  TstreamL1  = HitL1P * TwordL1r;
  TnstreamL1 = HitL1NP * TwordL1rNP;
  TtotalL1   = TstreamL1 + TnstreamL1;

  DEBUG_FPRINTF( stderr, "\tL1: regLoopR: %g, regLoopW: %g, Reads L1 Total: %g - Write L1 Total: %g, Kread: %d\n",
         regLoopR, regLoopW, regLoopR * I * J * K * NB, regLoopW * I * J * K * NB, Kread);
  DEBUG_FPRINTF( stderr, "\tL1-MissP: %.0lf (clines), L1-MissNP: %.0lf (clines), L1-HitP: %.0lf (words), L1-HitNP: %.0lf (words)\n",
         MissL1P, MissL1NP, HitL1P, HitL1NP);
  DEBUG_FPRINTF( stderr, "\tL1: nplanesL1: %g, nplanesL1NP: %g, nplanesL1P: %g\n", nplanesL1, nplanesL1NP, nplanesL1P);
  DEBUG_FPRINTF( stderr, "\tL1: TtotalL1: %g, TstreamL1: %g, TnstreamL1: %g, TwordL1r: %g\n", TtotalL1, TstreamL1, TnstreamL1, TwordL1r);


  /*****************************************/
  /* L1 - L2: Cost transfer from L2 to CPU */
  /*****************************************/
  DEBUG_FPRINTF( stderr, "*** L1 - L2 ***\n" );
  TclineL2r = cacheline / bwL2r;
  //TclineL2w = cacheline / bwL2w;
  TclineL2rNP = cacheline / bwL2rNP;
  //TclineL2wNP = cacheline / bwL2wNP;


  /* Compute Read Misses - L2 */
  nplanesL2 = read_misses( "L2", I, II, J, K, length, sizeL2, assocL2, INTERP, &Kread, &writePlanes, &nplanesL2NCI );

  volumeRead    = II * JJ * Kread * NB;                             // Read volume in points
  volumeReadCL  = ceil((II * JJ * Kread) / elemCLine) * NB;         // Read volume in Cache Lines


  /* Compute prefetching efficiency for L2 */
  /* For architectures with Prefetching system - Store reads may be prefetched as well - WRITEBACK */
  if ((prefL2) && (WRITEBACK)) {
    /* Add writePlanes because Writeback prefetch */
    //y0 = prefEL2[MIN(((int)nplanesL2)+((int)writePlanes), nPrefEffL2)-1];
    //y1 = prefEL2[MIN(((int)nplanesL2)+((int)writePlanes)+1, nPrefEffL2)-1];
    y0 = prefEL2[MIN(((int)nplanesL2)+((int)writePlanes), nPrefEffL2)-1];
    y1 = prefEL2[MIN(((int)nplanesL2)+((int)writePlanes)+1, nPrefEffL2)-1];
    //y0 = prefEL2[MIN(((int)nplanesL2NCI)+((int)writePlanes), nPrefEffL2)-1];
    //y1 = prefEL2[MIN(((int)nplanesL2NCI)+((int)writePlanes)+1, nPrefEffL2)-1];
    //y0 = prefEL2[MIN(((int)ceil(nplanesL2NCI))+((int)writePlanes), nPrefEffL2)-1];
    //y1 = prefEL2[MIN(((int)ceil(nplanesL2NCI))+((int)writePlanes)+1, nPrefEffL2)-1];
    //y0 = prefEL2[MIN(((int)nplanesL2), nPrefEffL2)-1];
    //y1 = prefEL2[MIN(((int)nplanesL2)+1, nPrefEffL2)-1];
  }
  else {
    y0 = prefEL2[MIN(((int)nplanesL2),   nPrefEffL2)-1];
    y1 = prefEL2[MIN(((int)nplanesL2)+1, nPrefEffL2)-1];
  }
  pEL2 = interpolate( INTERP, nplanesL2+((int)writePlanes), ((int)nplanesL2)+((int)writePlanes), ((int)nplanesL2)+((int)writePlanes)+1, y0, y1 );
  //pEL2 = interpolate( INTERP, nplanesL2NCI+((int)writePlanes), ((int)nplanesL2NCI)+((int)writePlanes), ((int)nplanesL2NCI)+((int)writePlanes)+1, y0, y1 );


  /* Blocking effect in Prefetched planes */
  if (NB > 1) {
    if (II > TP*elemCLine) { // Trigger Prefetching activation
      //I  += LAP*elemCLine; // Look-ahead Prefetching size
      //II += LAP*elemCLine; // Look-ahead Prefetching size
      //prefBlocking = ((TP+LAP)*elemCLine)/(II+((TP+LAP)*elemCLine));
      prefBlocking = 1.0;
      //if (nplanesL2P != 0.0) nplanesL2P += nplanesL2P * prefBlocking;
      //volumeRead   = (II+((TP+LAP)*elemCLine)) * JJ * Kread * NB;
      //volumeReadCL  = ceil(((II+((TP+LAP)*elemCLine)) * JJ * Kread) / elemCLine) * NB;         // Read volume in Cache Lines
    }
    else pEL2 = prefEL2[nPrefEffL2-1]; // Prefetching is not triggered
  }
  else prefBlocking = 1.0;
  //DEBUG_FPRINTF( stderr, "prefBlocking: %f\n", prefBlocking);
  //pEL2 *= prefBlocking;


#if 1
  nplanesL2P  = nplanesL2 * pEL2;         // Planes to read WITH prefetching
  nplanesL2NP = nplanesL2 * (1.0 - pEL2); // Planes to read WITHOUT prefetching
#else
  nplanesL2NP = MAX(nplanesL2 - prefL2[MIN(nplanesL2,nPrefEffL2)-1], 0.0); // Planes to read WITHOUT prefetching
  nplanesL2P  = nplanesL2 - nplanesL2NP;                                 // Planes to read WITH prefetching
#endif


  /* If NB and Trigger prefetching activated, add more CL due to Triggering and Look-ahead effect */
  if ((NB > 1) && (II > TP*elemCLine)) {
    MissL2P  = nplanesL2P * volumeReadCL + (nplanesL2P * LAP * JJ * Kread * NB);
    MissL2NP = nplanesL2NP * volumeReadCL + (nplanesL2 * TP * JJ * Kread * NB);
  }
  else {
    MissL2P  = nplanesL2P * volumeReadCL;
    MissL2NP = nplanesL2NP * volumeReadCL;
  }
#if 0
  HitL2P   = MAX(MissL1P - MissL2P, 0);
  HitL2NP  = MAX(MissL1NP - MissL2NP, 0);
#else
  HitL2P   = MAX(nplanesL1 - nplanesL2, 0) * volumeReadCL * pEL1;
  HitL2NP  = MAX(nplanesL1 - nplanesL2, 0) * volumeReadCL * (1 - pEL1);
#endif

  TstreamL2  = HitL2P * TclineL2r;
  TnstreamL2 = HitL2NP * TclineL2rNP;
  TtotalL2   = TstreamL2 + TnstreamL2;

  DEBUG_FPRINTF( stderr, "\tL2-MissP: %.0lf, L2-MissNP: %.0lf, L2-HitP: %.0lf L2-HitNP: %.0lf\n", MissL2P, MissL2NP, HitL2P, HitL2NP);
  DEBUG_FPRINTF( stderr, "\tL2: nplanesL2: %g, nplanesL2NP: %g, nplanesL2P: %g, prefEL2: %g\n", nplanesL2, nplanesL2NP, nplanesL2P, pEL2);
  DEBUG_FPRINTF( stderr, "\tL2: TtotalL2: %g, TstreamL2: %g, TnstreamL2: %g, TclineL2r: %g\n", TtotalL2, TstreamL2, TnstreamL2, TclineL2r);


#if defined(L3PRESENT)
  /* L2 - L3: Cost transfer from L3 to CPU */
  TclineL3r = cacheline / bwL3r;
  //TclineL3w = cacheline / bwL3w;
  TclineL3rNP = cacheline / bwL3rNP;
  //TclineL3wNP = cacheline / bwL3wNP;


  /* Compute Read Misses - L3 */
  nplanesL3 = read_misses( "L3", I, II, J, K, length, sizeL3, assocL3, INTERP, &Kread, &writePlanes, &nplanesL3NCI);

  volumeRead    = II * JJ * Kread * NB;                             // Read volume in points
  volumeReadCL  = ceil((II * JJ * Kread) / elemCLine) * NB;         // Read volume in Cache Lines

  nplanesL3NP = MAX(nplanesL3 - prefL3, 0.0); // Planes to read WITHOUT prefetching
  nplanesL3P  = nplanesL3 - nplanesL3NP;      // Planes to read WITH prefetching

  MissL3P  = nplanesL3P * volumeReadCL;
  MissL3NP = nplanesL3NP * volumeReadCL;
#if 0
  HitL3P   = MAX(MissL2P - MissL3P, 0);
  HitL3NP  = MAX(MissL2NP - MissL3NP,0);
#else
  HitL3P   = MAX(nplanesL2 - nplanesL3, 0) * volumeReadCL * pEL2;
  HitL3NP  = MAX(nplanesL2 - nplanesL3, 0) * volumeReadCL * (1 - pEL2);
#endif

  TstreamL3  = HitL3P * TclineL3r;
  TnstreamL3 = HitL3NP * TclineL3rNP;
  TtotalL3   = TstreamL3 + TnstreamL3;

  DEBUG_FPRINTF( stderr, "***********************************HITL3: %lf %lf %lf\n", HitL3P, MissL2P, MissL3P);
  DEBUG_FPRINTF( stderr, "L3-MissP: %.0lf, L3-MissNP: %.0lf, L3-HitP: %.0lf L3-HitNP: %.0lf\n", MissL3P, MissL3NP, HitL3P, HitL3NP);
  DEBUG_FPRINTF( stderr, "L3: nplanesL3: %g, nplanesL3NP: %g, nplanesL3P: %g\n", nplanesL3, nplanesL3NP, nplanesL3P);
//  DEBUG_FPRINTF( stderr, "L3: nplanesL3: %g, ncolumnsL3: %g, planesInter: %g, planesReuse: %g\n", nplanesL3, ncolumnsL3, planesInter, planesReuse);


  /* L3 - RAM: Cost transfer from MEM to L3 those not streamed and missed */
  TclineRAMr = cacheline / bwRAMr;
  TclineRAMw = cacheline / bwRAMw;
  TclineRAMrNP = cacheline / bwRAMrNP;
  //TclineRAMwNP = cacheline / bwRAMwNP;

  TstreamRAM  = MissL3P * TclineRAMr;
  TnstreamRAM = MissL3NP * TclineRAMrNP;
  TtotalRAM   = TstreamRAM + TnstreamRAM;


  /************** WRITE ***************************/
  /* Write data back into memory I * J * K volume */
#if defined(WRITEBACK)
  //Twrite = regLoopW * I * J * K * NB * TwordL1w + writePlanes * volumeWriteCL * TclineRAMwNP;
  nplanesWP  = writePlanes * pEL3;
  nplanesWNP = writePlanes * (1.0 - pEL3);
  MissWP     = nplanesWP  * volumeWriteCL;
  MissWNP    = nplanesWNP * volumeWriteCL;

  TwriteP  = MissWP  * TclineRAMw;
  TwriteNP = MissWNP * TclineRAMwNP;
  Twrite   = TwriteP + TwriteNP;
  DEBUG_FPRINTF( stderr, "\tWRITEBACK - Twrite: %g, volumeWriteCL: %g, TwriteP: %g, TwriteNP: %g\n",
               Twrite, volumeWriteCL, TwriteP, TwriteNP);
#else
  Twrite = regLoopW * I * J * K * NB * TwordL1w + writePlanes * volumeWriteCL * TclineRAMw;
  DEBUG_FPRINTF( stderr, "\tWRITETHROUGH - Twrite: %g, volumeWriteCL: %g, TclineRAMw: %g, WriteOutput: %g\n",
               Twrite, volumeWrite, TclineRAMw, volumeWriteCL * TclineRAMw);
#endif


  Ttotal      = Twrite + /*Tcold +*/ TtotalL1 + TtotalL2 + TtotalL3 + TtotalRAM;
  DEBUG_FPRINTF( stderr, "RAM-HitP: %g RAM-HitNP: %g\n", MissL3P, MissL3NP);

  //DEBUG_FPRINTF( stderr, "%d %d %d %d %d %d %d %d %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %0.lf\n", nx, ny, nz, tx, ty, tz, timesteps, length, MissL1P, MissL1NP, HitL1P, HitL1NP, MissL2P, MissL2NP, HitL2P, HitL2NP, MissL3P, MissL3NP, HitL3P, HitL3NP );
#else

  /* L2 - RAM: Cost transfer from MEM to L2 those not streamed and missed */
  DEBUG_FPRINTF( stderr, "*** L2 - MEM ***\n" );
  TclineRAMr = cacheline / bwRAMr;
  TclineRAMw = cacheline / bwRAMw;
  TclineRAMrNP = cacheline / bwRAMrNP;
  TclineRAMwNP = cacheline / bwRAMwNP;


  TstreamRAM = MissL2P * TclineRAMr;
  TnstreamRAM = MissL2NP * TclineRAMrNP;
  TtotalRAM  = TstreamRAM + TnstreamRAM;
  DEBUG_FPRINTF( stderr, "\tMissToRAM: MissL2P: %g, MissL2NP: %g\n", MissL2P, MissL2NP);
  DEBUG_FPRINTF( stderr, "\tRAM: TtotalRAM: %g, TstreamRAM: %g, TnstreamRAM: %g\n", TtotalRAM, TstreamRAM, TnstreamRAM);


  /* Bring cold data from MEM to CPU */
  //Tcold      = ReadPlanes * ceil((II * JJ) / elemCLine) * TclineRAMr;
  //DEBUG_FPRINTF( stderr, "Tcold: %g\n", Tcold/TclineRAMr);


  /************** WRITE ***************************/
  /* Write data back into memory I * J * K volume */
#if defined(WRITEBACK)
  //Twrite = regLoopW * I * J * K * NB * TwordL1w + writePlanes * volumeWriteCL * TclineRAMwNP;
  nplanesWP  = writePlanes * pEL2;
  nplanesWNP = writePlanes * (1.0 - pEL2);
  MissWP     = nplanesWP  * volumeWriteCL;
  MissWNP    = nplanesWNP * volumeWriteCL;

  TwriteP  = MissWP  * TclineRAMw;
  TwriteNP = MissWNP * TclineRAMwNP;
  Twrite   = TwriteP + TwriteNP;
  DEBUG_FPRINTF( stderr, "\tWRITEBACK - writePlanes: %g, Twrite: %g, volumeWriteCL: %g, TwriteP: %g, TwriteNP: %g\n",
               writePlanes, Twrite, volumeWriteCL, TwriteP, TwriteNP);
#else
  Twrite = regLoopW * I * J * K * NB * TwordL1w + writePlanes * volumeWriteCL * TclineRAMw;
  DEBUG_FPRINTF( stderr, "\tWRITETHROUGH - Twrite: %g, volumeWriteCL: %g, TclineRAMw: %g, WriteOutput: %g\n",
               Twrite, volumeWrite, TclineRAMw, volumeWriteCL * TclineRAMw);
#endif


  Ttotal = Twrite + /*Tcold +*/ TtotalL1 + TtotalL2 + TtotalRAM;
  DEBUG_FPRINTF( stderr, "*** Ttotal: %g, Twrite: %g, TtotalL1: %g, TtotalL2: %g, TtotalRAM: %g\n",
                Ttotal, Twrite, TtotalL1, TtotalL2, TtotalRAM);
  //DEBUG_FPRINTF( stderr, "%d %d %d %d %d %d %d %d %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf\n",
  //        nx, ny, nz, tx, ty, tz, timesteps, length,
  //        MissL1P, MissL1NP, HitL1P, MissL2P, MissL2NP, HitL2P, HitL2NP );

#endif

  DEBUG_FPRINTF( stderr, "IDEAL (CASE 3): ######################################\n" );


#if 0
/* TODO: Timesteps is printed but not used to update values - Run with timesteps 1 */
  //printf("#Title NX NY NZ TX TY TZ Timesteps Lenght"
  //       " nplanesL1P nplanesL1NP MissL1P MissL1NP HitL1P HitL1NP prefEL1" 
  //       " nplanesL2P nplanesL2NP MissL2P MissL2NP HitL2P HitL2NP prefEL2"
  //       " nplanesL3P nplanesL3NP MissL3P MissL3NP HitL3P HitL3NP prefEL3"
  //       " nplanesWP  nplanesWNP  MissWP  MissWNP"  
  //       " TtotalL1 TtotalL2 TtotalL3 TtotalRAM Twrite Ttotal\n");
  fprintf( stdout, "\"%dx%d\" %d %d %d %d %d %d %d %d"
           " %.10g %.10g %.10g %.10g %.10g %.10g %.10g"
           " %.10g %.10g %.10g %.10g %.10g %.10g %.10g"
           " %.10g %.10g %.10g %.10g %.10g %.10g %.10g"
           " %.10g %.10g %.10g %.10g"
           " %.10g %.10g %.10g %.10g %.10g %.10g\n",
           tx, ty, nx, ny, nz, tx, ty, tz, timesteps, length,
           nplanesL1P, nplanesL1NP, MissL1P, MissL1NP, HitL1P, HitL1NP, pEL1,
           nplanesL2P, nplanesL2NP, MissL2P, MissL2NP, HitL2P, HitL2NP, pEL2,
           nplanesL3P, nplanesL3NP, MissL3P, MissL3NP, HitL3P, HitL3NP, pEL3,
           nplanesWP,  nplanesWNP,  MissWP,  MissWNP,
           TtotalL1, TtotalL2, TtotalL3, TtotalRAM, Twrite, Ttotal );
#endif

  /* Return execution time */
  return Ttotal;
}


/* 
 * Command-line parameters
 *
 * argv[1] -> nx
 * argv[2] -> ny
 * argv[3] -> nz
 * argv[4] -> tx
 * argv[5] -> ty
 * argv[6] -> tz
 * argv[7] -> timesteps
 * argv[8] -> length
 */
int main( int argc, char **argv )
{
  /* input */
  int nx, ny, nz;
  int tx, ty, tz;
  int timesteps, length;
  /* output */
  double to, tp, ti;

  /* checking command line parameters */
  if( argc < 9 ) 
  {
    fprintf( stderr, "Usage:\n" );
    fprintf( stderr, "%s nx ny nz tx ty tz timesteps length\n", argv[0] );
    return EXIT_FAILURE;
  }

#ifdef HAVE_PAPI
  PAPI_library_init(PAPI_VER_CURRENT);
#endif

  /* reading command line parameters */
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  tx = atoi(argv[4]);
  ty = atoi(argv[5]);
  tz = atoi(argv[6]);
  timesteps = atoi(argv[7]);
  length = atoi(argv[8]);

  /* 
   * Assuming computation is fully overlap with memory 
   * transfers, basically memory bound assumption. The 
   * following functions compute what will be memory movement time. 
   * */
  /* optimistic implies: ? */
  to = transfer_optimistic( nx, ny, nz, length, timesteps );
  /* pesimistic implies: ? */
  tp = transfer_pesimistic( nx, ny, nz, length, timesteps );
  /* ideal stands for: ? */
  ti = transfer_ideal( nx, ny, nz, tx, ty, tz, length, timesteps );

  __FPRINTF( stderr, "nx*ny*nz \ttx*ty*tz \ttimesteps \tlength"
                 " \toptimistic \tpesimistic \tideal \n" );
  __FPRINTF( stderr, "%d %d %d \t%d %d %d \t%d \t\t%d \t%g \t%g \t%g\n",
                 nx, ny, nz, tx, ty, tz, timesteps, length, to, tp, ti );

  return EXIT_SUCCESS;
}
