
/* Number of Threads per Core */
#ifndef THRCORE
# define THRCORE 1
#endif

/* Common definitions */
#define KB (1024)
#define MB (1024*1024)
#define GB (1024*1024*1024)

/* Architectureal definitions */

#if NEHALEM
# define bwL1r      (49.4*1024*1024*1024) /* (32.6*1024*1024*1024) */
# define bwL2r      (29.4*1024*1024*1024)
# define bwL3r      (21.1*1024*1024*1024)
# define bwRAMr     (8.2*1024*1024*1024)

# define bwL1w      (49.4*1024*1024*1024)
# define bwL2w      (30.5*1024*1024*1024)
# define bwL3w      (17.6*1024*1024*1024)
# define bwRAMw     (7.9*1024*1024*1024)

# define bwL1rNP    (23.4*1024*1024*1024)
# define bwL2rNP    (12.5*1024*1024*1024)
# define bwL3rNP    (8.1*1024*1024*1024)
# define bwRAMrNP   (3.5*1024*1024*1024)

# define bwL1wNP    (24.6*1024*1024*1024)
# define bwL2wNP    (15.3*1024*1024*1024)
# define bwL3wNP    (9.5*1024*1024*1024)
# define bwRAMwNP   (3.9*1024*1024*1024)

# define cacheline  64.0
# define assocL1    8.0
# define assocL2    8.0
# define assocL3    16.0
# define sizeL1     (32.0*1024)
# define sizeL2     (256.0*1024)
# define sizeL3     (8.0*1024*1024)  /* shared 4 cores */
# define prefL1     2.0
# define prefL2     2.0
# define prefL3     0.0
# define regFPR     12.0 /* http://en.wikipedia.org/wiki/AMD64 */
# define regGPR     10.0
# define L3PRESENT  1
# define WRITEBACK  1

#elif SANDY
# define bwL1r      (49.4*1024*1024*1024) /* (32.6*1024*1024*1024) */
# define bwL2r      (29.4*1024*1024*1024)
# define bwL3r      (21.1*1024*1024*1024)
# define bwRAMr     (8.2*1024*1024*1024)

# define bwL1w      (49.4*1024*1024*1024)
# define bwL2w      (30.5*1024*1024*1024)
# define bwL3w      (17.6*1024*1024*1024)
# define bwRAMw     (7.9*1024*1024*1024)

# define bwL1rNP    (23.4*1024*1024*1024)
# define bwL2rNP    (12.5*1024*1024*1024)
# define bwL3rNP    (8.1*1024*1024*1024)
# define bwRAMrNP   (3.5*1024*1024*1024)

# define bwL1wNP    (24.6*1024*1024*1024)
# define bwL2wNP    (15.3*1024*1024*1024)
# define bwL3wNP    (9.5*1024*1024*1024)
# define bwRAMwNP   (3.9*1024*1024*1024)

# define cacheline  64.0
# define assocL1    8.0
# define assocL2    8.0
# define assocL3    16.0
# define sizeL1     (32.0*1024)
# define sizeL2     (256.0*1024)
# define sizeL3     (8.0*1024*1024)  /* shared 4 cores */
# define prefL1     2.0
# define prefL2     2.0
# define prefL3     0.0
# define regFPR     12.0 /* http://en.wikipedia.org/wiki/AMD64 */
# define regGPR     10.0
# define L3PRESENT  1
# define WRITEBACK  1

# define nPrefEffL1 1   // Max index in prefEffL1
# define nPrefEffL2 13  // Max index in prefEffL2

# define prefEffL1  {0.0}
# define prefEffL2  {0.983, 0.983, 0.983, 0.983, 0.975, 0.973, 0.963, 0.952, 0.945, 0.929, 0.853, 0.210, 0.031}
# define VECTOR     0 // Vector code disabled
# define UNALIGNED  1 // Unaligned read/write?

# define TP         3 // Triggering Prefetching (in # of clines)
# define LAP        5 // Look-ahead Prefetching (in # of clines)

#elif IVY // My Acer S7
# define bwL1r      (30.0*1024*1024*1024)
# define bwL2r      (27.5*1024*1024*1024)
# define bwL3r      (24.9*1024*1024*1024)
# define bwRAMr     (14.4*1024*1024*1024)

# define bwL1w      (46.0*1024*1024*1024)
# define bwL2w      (29.8*1024*1024*1024)
# define bwL3w      (21.5*1024*1024*1024)
# define bwRAMw     (7.4*1024*1024*1024)

# define bwL1rNP    (20.6*1024*1024*1024)
# define bwL2rNP    (14.9*1024*1024*1024)
# define bwL3rNP    (8.5*1024*1024*1024) //???
# define bwRAMrNP   (3.7*1024*1024*1024)

# define bwL1wNP    (21.0*1024*1024*1024)
# define bwL2wNP    (18.9*1024*1024*1024)
# define bwL3wNP    (15.3*1024*1024*1024)
# define bwRAMwNP   (2.6*1024*1024*1024)

# define cacheline  64.0
# define assocL1    8.0
# define assocL2    8.0
# define assocL3    8.0
# define sizeL1     (32.0*1024)
# define sizeL2     (256.0*1024)
# define sizeL3     (4.0*1024*1024)
# define prefL1     16.0
# define prefL2     16.0
# define prefL3     16.0
# define regFPR     12.0 /* http://en.wikipedia.org/wiki/AMD64 */
# define regGPR     10.0
# define L3PRESENT  1
# define WRITEBACK  1

#elif KNC

/* BW performance must be set to the maximum Stream2
 * performance obtained when running with the number
 * of threads (cores) to perform the simulation */

#if THREADS == 1
# error "Not done!"

#elif THREADS == 4

/* 4 threads - 1 thread per core */
# if THRCORE == 1
# define bwL1r      (61.2*GB)
# define bwL2r      (11.3*GB)
# define bwRAMr     (8.55*GB)

# define bwL1w      (54.5*GB)
# define bwL2w      (10.3*GB)
# define bwRAMw     (4.6*GB) //(0.92*GB)

# define bwL1rNP    (54.9*GB) //(61.1*GB)
# define bwL2rNP    (11.0*GB) //(11.4*GB)
# define bwRAMrNP   (0.90*GB) //(0.91*GB)

# define bwL1wNP    (54.5*GB)
# define bwL2wNP    (11.3*GB)
# define bwRAMwNP   (0.91*GB) //(0.89*GB) //(0.92*GB)

/* 4 threads - 2 threads per core */
# elif THRCORE == 2
# define bwL1r      (46.5*GB)
# define bwL2r      (10.7*GB)
# define bwRAMr     (10.7*GB) // Prefetched data use previous L2 bw (6.33*GB)
                              
# define bwL1w      (55.9*GB)
# define bwL2w      (11.2*GB)
# define bwRAMw     (11.2*GB) // Prefetched data use previous L2 bw (3.85*GB)

# define bwL1rNP    (47.5*GB)
# define bwL2rNP    (10.8*GB)
# define bwRAMrNP   (0.90*GB)

# define bwL1wNP    (56.0*GB)
# define bwL2wNP    (11.2*GB)
# define bwRAMwNP   (0.91*GB)

/* 4 threads - 4 threads per core */
# elif THRCORE == 4
# define bwL1r      (28.0*GB)
# define bwL2r      (10.4*GB)
# define bwRAMr     (0.90*GB)
                              
# define bwL1w      (27.4*GB)
# define bwL2w      (10.6*GB)
# define bwRAMw     (0.91*GB)

# define bwL1rNP    (28.0*GB)
# define bwL2rNP    (10.4*GB)
# define bwRAMrNP   (0.90*GB)

# define bwL1wNP    (27.4*GB)
# define bwL2wNP    (10.6*GB)
# define bwRAMwNP   (0.91*GB)
# else
#   error "Not done!"
#endif
#elif THREADS == 16

# if THRCORE == 1
# define bwL1r      (50.0*GB)//(28.0*GB)
# define bwL2r      (39.0*GB) //(24.0*GB)
# define bwRAMr     (3.5*GB)  //(3.5*GB)

# define bwL1w      (290.0*GB)
# define bwL2w      (38.4*GB)
# define bwRAMw     (3.5*GB)

# define bwL1rNP    (50.0*GB)
# define bwL2rNP    (39.0*GB)
# define bwRAMrNP   (3.5*GB)

# define bwL1wNP    (290.0*GB)
# define bwL2wNP    (38.4*GB)
# define bwRAMwNP   (3.5*GB)
# elif THRCORE == 2
# define bwL1r      (28.0*GB) //(80.0*GB)
# define bwL2r      (24.0*GB) //(42.1*GB)
# define bwRAMr     (3.5*GB)  //(3.4*GB)
                               
# define bwL1w      (290.0*GB)//(235.0*GB)
# define bwL2w      (38.4*GB) //(44.0*GB)
# define bwRAMw     (3.5*GB)  //(3.5*GB)
                               
# define bwL1rNP    (28.0*GB) //(80.0*GB)
# define bwL2rNP    (24.0*GB) //(42.1*GB)
# define bwRAMrNP   (3.5*GB)  //(3.4*GB)
                               
# define bwL1wNP    (290.0*GB)//(235.0*GB)
# define bwL2wNP    (38.4*GB) //(44.0*GB)
# define bwRAMwNP   (3.5*GB)  //(3.5*GB)
# elif THRCORE == 4
# define bwL1r      (28.0*GB) //(80.0*GB)
# define bwL2r      (24.0*GB) //(41.3*GB)
# define bwRAMr     (3.5*GB)  //(3.5*GB)
                                         
# define bwL1w      (290.0*GB)//(105.0*GB)
# define bwL2w      (38.4*GB) //(42.6*GB)
# define bwRAMw     (3.5*GB)  //(3.5*GB)
                              
# define bwL1rNP    (28.0*GB) //(80.0*GB)
# define bwL2rNP    (24.0*GB) //(41.3*GB)
# define bwRAMrNP   (3.5*GB)  //(3.5*GB)
                              
# define bwL1wNP    (290.0*GB)//(105.0*GB)
# define bwL2wNP    (38.4*GB) //(42.6*GB)
# define bwRAMwNP   (3.5*GB)  //(3.5*GB)
# else
#   error "Not done!"
# endif

#elif THREADS == 32
# define bwL1r      (215.0*GB)
# define bwL2r      (74.0*GB)
# define bwRAMr     (7.0*GB)
                                  
# define bwL1w      (560.0*GB)
# define bwL2w      (78.4*GB)
# define bwRAMw     (6.98*GB)

# define bwL1rNP    (215.0*GB)
# define bwL2rNP    (74.0*GB)
# define bwRAMrNP   (7.0*GB)

# define bwL1wNP    (560.0*GB)
# define bwL2wNP    (78.4*GB)
# define bwRAMwNP   (6.98*GB)
#elif THREADS == 61

/* 61 threads - 1 thread per core */
# define bwL1r      (115.0*GB)
# define bwL2r      (89.0*GB)
# define bwRAMr     (12.6*GB)

# define bwL1w      (700.0*GB)
# define bwL2w      (100.4*GB)
# define bwRAMw     (11.4*GB)

# define bwL1rNP    (115.0*GB)
# define bwL2rNP    (89.0*GB)
# define bwRAMrNP   (12.6*GB)

# define bwL1wNP    (700.0*GB)
# define bwL2wNP    (100.4*GB)
# define bwRAMwNP   (11.4*GB)
#else
# error "Number of threads not done!"
#endif


#if THRCORE == 1

/* Data obtained from: prefetchers.noprefetch.mic.1thrxcore.perf3 */
# define nPrefEffL1 1   // Max index in prefEffL1
# define nPrefEffL2 22  // Max index in prefEffL2

// For a higher number of prefetching streams we use the last coefficient
//                   1stre  2stre  3stre  4stre  5stre  6stre  7stre  8stre  9stre  10stre 11stre 12str  13str  14str  15str  16str, 17str, 18str, 19str, 20str, 21str, 22str
# define prefEffL1  {0.0}
# define prefEffL2  {0.983, 0.983, 0.983, 0.983, 0.982, 0.983, 0.983, 0.982, 0.981, 0.980, 0.981, 0.976, 0.976, 0.968, 0.979, 0.965, 0.961, 0.942, 0.073, 0.075, 0.053, 0.031}

#elif THRCORE == 2

/* Data obtained from: prefetchers.noprefetch.mic.2thrxcore.perf3 */
# define nPrefEffL1 1   // Max index in prefEffL1
# define nPrefEffL2 13  // Max index in prefEffL2

# define prefEffL1  {0.0}
# define prefEffL2  {0.983, 0.983, 0.983, 0.983, 0.975, 0.973, 0.963, 0.952, 0.945, 0.929, 0.853, 0.210, 0.031}

#elif THRCORE == 3

# error "Not yet done!"

#elif THRCORE == 4

/* Data obtained from: prefetchers.noprefetch.mic.4thrxcore.perf3 */
# define nPrefEffL1 1   // Max index in prefEffL1
# define nPrefEffL2 7   // Max index in prefEffL2

# define prefEffL1  {0.0}
# define prefEffL2  {0.981, 0.981, 0.970, 0.953, 0.930, 0.900, 0.031}

#else
# error "THRCORE not supported for MIC platform"
#endif


/* Evaluate TLB misses cost!!!! */

/* Cache sizes only must be modified when
 * several threads are running on the same
 * core, otherwise they keep the whole sizes. */
# define cacheline  64.0
# define assocL1    (8.0/THRCORE)
# define assocL2    (8.0/THRCORE)
# define sizeL1     ((32.0*KB)/(THRCORE))
# define sizeL2     ((512.0*KB)/(THRCORE))
# define prefL1     0.0
# define prefL2     16.0 // 16 HW prefetchers per core (shared among threads)
# define regFPR     26.0 //32.0 //12.0 /* 32 Vector registers Intel Xeon Phi FIGURE 8.7 */
# define regGPR     8.0 //10.0 /* 16 GPR but some of them are used for other things */
# define WRITEBACK  1
# define VECTOR     1 // Vector code enabled
# define VECSIZE    8 // Vector code enabled
# define UNALIGNED  1 // Unaligned read/write?

# define TP 3         // Triggering Prefetching (in # of clines)
# define LAP 5        // Look-ahead Prefetching (in # of clines)

#elif OPTERON

# define bwL1r      (29.1*1024*1024*1024)
# define bwL2r      (13.9*1024*1024*1024)
# define bwL3r      (7.6*1024*1024*1024)
# define bwRAMr     (3.3*1024*1024*1024)

# define bwL1w      (29.9*1024*1024*1024)
# define bwL2w      (12.1*1024*1024*1024)
# define bwL3w      (4.9*1024*1024*1024)
# define bwRAMw     (4.9*1024*1024*1024)

# define bwL1rNP    (14.6*1024*1024*1024)
# define bwL2rNP    (4.1*1024*1024*1024)
# define bwL3rNP    (2.2*1024*1024*1024)
# define bwRAMrNP   (1.3*1024*1024*1024)

# define bwL1wNP    (8.6*1024*1024*1024)
# define bwL2wNP    (4.1*1024*1024*1024)
# define bwL3wNP    (2.1*1024*1024*1024)
# define bwRAMwNP   (0.8*1024*1024*1024)

# define cacheline  64.0
# define assocL1    2.0
# define assocL2    16.0
# define assocL3    32.0
# define sizeL1     (64.0*1024)
# define sizeL2     (512.0*1024)
# define sizeL3     (2.0*1024*1024)  /* shared */
# define prefL1     2.0 /*2   // 8 prefetchers for all the chip */
# define prefL2     0.0
# define prefL3     0.0
# define regFPR     12.0 /* http://en.wikipedia.org/wiki/AMD64 */
# define regGPR     10.0
# define L3PRESENT  1
# define WRITEBACK  1

#elif SHANGHAI
# define bwL1r      (34.1*1024*1024*1024)
# define bwL2r      (15.9*1024*1024*1024)
# define bwL3r      (8.9*1024*1024*1024)
# define bwRAMr     (2.8*1024*1024*1024)

# define bwL1w      (42.2*1024*1024*1024)
# define bwL2w      (14.6*1024*1024*1024)
# define bwL3w      (5.0*1024*1024*1024)
# define bwRAMw     (5.0*1024*1024*1024)

# define bwL1rNP    (16.6*1024*1024*1024)
# define bwL2rNP    (5.1*1024*1024*1024)
# define bwL3rNP    (2.7*1024*1024*1024)
# define bwRAMrNP   (1.2*1024*1024*1024)

# define bwL1wNP    (9.9*1024*1024*1024)
# define bwL2wNP    (4.9*1024*1024*1024)
# define bwL3wNP    (2.9*1024*1024*1024)
# define bwRAMwNP   (1.0*1024*1024*1024)

# define cacheline  64.0 /* Architecture of the AMD Quad Core CPUs Brian
# Waldecker, Ph.D. | April 13, 2009 */
# define assocL1    2.0
# define assocL2    16.0
# define assocL3    48.0
# define sizeL1     (64.0*1024)
# define sizeL2     (512.0*1024)
# define sizeL3     (6.0*1024*1024)  /* shared */
# define prefL1     2.0 /*2   // 8 prefetchers for all the chip */
# define prefL2     0.0
# define prefL3     0.0
# define regFPR     12.0 /* http://en.wikipedia.org/wiki/AMD64 */
# define regGPR     10.0
# define L3PRESENT  1
# define WRITEBACK  1

#elif POWER6
# define bwL1r      (51.8*1024*1024*1024)
# define bwL2r      (35.4*1024*1024*1024)
# define bwL3r      (16.6*1024*1024*1024)
# define bwRAMr     (8.2*1024*1024*1024)

# define bwL1w      (20.7*1024*1024*1024)
# define bwL2w      (20.7*1024*1024*1024)
# define bwL3w      (12.3*1024*1024*1024)
# define bwRAMw     (4.3*1024*1024*1024)

# define bwL1rNP    (22.0*1024*1024*1024)
# define bwL2rNP    (19.5*1024*1024*1024)
# define bwL3rNP    (15.1*1024*1024*1024)
# define bwRAMrNP   (7.3*1024*1024*1024)

# define bwL1wNP    (14.1*1024*1024*1024)
# define bwL2wNP    (13.9*1024*1024*1024)
# define bwL3wNP    (11.2*1024*1024*1024)
# define bwRAMwNP   (3.9*1024*1024*1024)

# define cacheline  128.0
# define assocL1    8.0
# define assocL2    8.0
# define assocL3    16.0
# define sizeL1     (64.0*1024)
# define sizeL2     (4.0*1024*1024)
# define sizeL3     (32.0*1024*1024)  /* shared */
# define prefL1     8.0
# define prefL2     8.0
# define prefL3     0.0
# define regFPR     26.0 /* http://en.wikipedia.org/wiki/PowerPC */
# define regGPR     26.0
# define L3PRESENT  1
# define WRITETHROUGH 1

#elif BGP
# define bwL1r      (6.2*1024*1024*1024)
# define bwL2r      (2.1*1024*1024*1024)
# define bwL3r      (2.0*1024*1024*1024)
# define bwRAMr     (2.0*1024*1024*1024)

# define bwL1w      (3.3*1024*1024*1024)
# define bwL2w      (3.3*1024*1024*1024)
# define bwL3w      (3.3*1024*1024*1024)
# define bwRAMw     (3.3*1024*1024*1024)

# define bwL1rNP    (3.2*1024*1024*1024)
# define bwL2rNP    (1.3*1024*1024*1024)
# define bwL3rNP    (0.6*1024*1024*1024)
# define bwRAMrNP   (0.6*1024*1024*1024)

# define bwL1wNP    (3.3*1024*1024*1024)
# define bwL2wNP    (3.3*1024*1024*1024)
# define bwL3wNP    (3.3*1024*1024*1024)
# define bwRAMwNP   (1.9*1024*1024*1024)

# define cacheline  128.0
# define cachelineL1 32.0  /* What's up with this??? */
# define cachelineL2 128.0
# define cachelineL3 128.0
# define assocL1    64.0
# define assocL2    1.0          /* Fully associative */
# define assocL3    8.0
# define sizeL1     (32.0*1024)  /* Write-through */
# define sizeL2     (1920)
# define sizeL3     (8.0*1024*1024)
# define prefL1     0.0
# define prefL2     7.0
# define prefL3     0.0
# define regFPR     26.0 /* http://en.wikipedia.org/wiki/PowerPC */
# define regGPR     26.0
# define L3PRESENT  1
# define WRITETHROUGH 1

#elif PPC970
# define bwL1r      (23.7*1024*1024*1024)
# define bwL1w      (6.1*1024*1024*1024)
# define bwL2r      (13.0*1024*1024*1024)
# define bwL2w      (2.8*1024*1024*1024)
# define bwL3r      (0.0*1024*1024*1024)
# define bwL3w      (0.0*1024*1024*1024)
# define bwRAMr     (4.9*1024*1024*1024)
# define bwRAMw     (1.5*1024*1024*1024)
# define cacheline  128.0
# define assocL1    2.0
# define assocL2    8.0
# define assocL3    0.0
# define sizeL1     (32.0*1024)
# define sizeL2     (1.0*1024*1024)
# define sizeL3     (0.0*1024*1024)  /* shared */
# define prefL1     0.0
# define prefL2     0.0  /* 8.0 */
# define prefL3     0.0
# define regFPR     26.0 /* http://en.wikipedia.org/wiki/PowerPC */
# define regGPR     26.0
#endif
