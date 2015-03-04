

#define __PRINTF(fmt, ... )             printf(fmt, ## __VA_ARGS__)
#define __FPRINTF(file, fmt, ...)       fprintf(file, fmt, ## __VA_ARGS__)

#ifdef DEBUG
# define DEBUG_PRINTF(fmt, ...)         __PRINTF(fmt, ## __VA_ARGS__)
# define DEBUG_FPRINTF(file, fmt, ...)  __FPRINTF(file, fmt, ## __VA_ARGS__)
#else
# define DEBUG_PRINTF(fmt, ...)
# define DEBUG_FPRINTF(file, fmt, ...)
#endif


#define word 8 /* double */
#define Cops 0
#define dim 3.0 /* 3D problems? */

/* By default (0.5) -> Half of the cache
 * Otherwise columns central plane / columns in total */
#define COLUMN_RATIO (double)((2 * length + 1)/(2 * length * (dim - 1) + 1))
/* Add 1 because writing (write-back policy) */
//#define COLUMN_RATIO (double)((2 * length + 1)/((2 * length * (dim - 1) + 1) + 1))
//#define COLUMN_RATIO (double)((2 * length + 1)/(2 * length * dim + 1))

#define elemCLine (cacheline / word)

/* Not used in ideal */
#define bwL2  (14.0*1024*1024*1024)
#define bwRAM  (4.8*1024*1024*1024)
#define wordL2 ((6*1024*1024)/word)
#define REP 10
#define L2SIZE (6*1024*1024)
/*******************************/


#define MAX(x,y) (x > y ? x : y)
#define MIN(x,y) (x < y ? x : y)

/* Interpolation methods */
#define NO_INTERP       0
#define LIN_INTERP      1
#define EXP_INTERP      2
#define LOG_INTERP      3

/* By default no interpolation method */
#ifndef INTERP
# define INTERP NO_INTERP
#endif

