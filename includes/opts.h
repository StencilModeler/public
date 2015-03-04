/*
 * Plane and words read/written definitions
 * for stencils schemes and optimizations.
 */

#ifdef SEMI

/*
 *  o   -> Point to read
 *  O   -> Point to read/write
 *  R   -> Read
 *  W   -> Write
 *  t   -> Prev. timestep
 *  t+1 -> Next timestep
 *
 *         . R^t O
 *         .   /  W^t+1
 *       R^t o
 *         ./      R/W^t+1
 * R/W^t+1 O---o----O
 * R^t     |  R^t  R^t 
 *         .
 *         .
 *         .
 */

// Deprecated, should be removed at any time
# define ReadStencil  (double)((length * (dim - 1)) + 1)
# define ReadStencilL (double)(1)
# define ReadWordsL   (double)(1)
# define WriteStencil (double)(dim - 1)
// Deprecated, should be removed at any time

# if VECTOR == 1
                                /* Reads t              */   /* Reads t+1 */
#   define ReadWords   (double)(((length * (dim - 1)) + 1) + (dim - 1))
# else
                                /* Reads t                           */   /* Reads t+1 */
#   define ReadWords   (double)(((length * (dim - 1)) + 2 * length + 1) + (dim - 1))
# endif

                              /* length reads t + 2 reads t+1 */
# define ReadPlanes   (double)(length + 2)                   // Planes to read
# define WriteWords   (double)(3)                            // Words to write
/* 3 planes are written when ReadPlane*II is too big to fit in cache - no k-central plane reuse */
# define WritePlanes  (double)(2)                            // Planes to write: Low bound: 2 -> High bound: 3

// Classical stencil
#else

// Deprecated, should be removed at any time
# define ReadStencil  (double)((2 * length * (dim - 1)) + 1) // Cache lines to read
# define ReadStencilL (double)(2 * length + 1) // Planes to read
# define ReadWordsL   (double)((2 * length * (dim - 1)) + 1)
# define WriteStencil (double)(1)
// Deprecated, should be removed at any time

# if VECTOR == 1
#   define ReadWords  (double)((2 * length * (dim - 1)) + 1) // Central plane -> only read 1
# else
#   define ReadWords  (double)((2 * length * dim) + 1)       // Words to read
# endif

# define ReadPlanes   (double)(2 * length + 1)               // Planes to read
# define WriteWords   (double)(1)                            // Words to write
# define WritePlanes  (double)(1)                            // Planes to write

#endif

