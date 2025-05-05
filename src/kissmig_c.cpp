/* Michael Nobis & Dominik Landolt, May 2025 */

#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
Rcpp::NumericVector kissmig_c(const Rcpp::NumericVector &d, const Rcpp::NumericVector &s, const Rcpp::IntegerVector &dim, int it, double pe, double pc, uint8_t ty, int si, int nThreads, uint32_t randRingSize)
{
  int l;						// index of iteration
  uint32_t i,					// generic index of cell in layer
  t, 							// index of layer
  p,							// index of cell in layer
  p1,							// index of cell + 1 row
  p2,							// index of cell - 1 row
  flocIdx,						// index for FOC and LOC
  nrow,							// rows per layer
  ncol,							// columns per layer
  nlay,							// amount of layers
  n;							// amount of cells per layer
  uint8_t b = 1;				// flag for checking corner events
  const double* sPtr = &s[0];	// pointer to start of one suitability layer
  nrow = dim[0];
  ncol = dim[1];
  nlay = dim[2];
  n = nrow * ncol;
  Rcpp::NumericVector ans(n);				// resulting distribution, initializes with 0s
  uint8_t* arrT0 = R_Calloc(n, uint8_t); 	// temp distribution on previous time step
  uint8_t *arrT1 = R_Calloc(n, uint8_t);    // temp distribution on current time step
  uint32_t* occPos = R_Calloc(n, uint32_t); // stores positions of occurrences on current layer

  // R RNG
  double* randRing = R_Calloc(randRingSize, double); 	// store for random number used as ring buffer
  uint32_t idxRandomRing; 							    // index of ring buffer "randRing"
  // generate random number for ring buffer
  GetRNGstate();
  for (i = 0; i < randRingSize; i++)
  {
    randRing[i] = R::unif_rand();
  }


  uint32_t nIterationAndLayer = (nlay * it); // amount of total iterations including every layer
  uint32_t currIterationAndLayer = 0;        // tracking current iteration of total iterations

  uint32_t* idxRandomNumberStart = R_Calloc(nIterationAndLayer, uint32_t); // random index in random number per iterations in total
  // generate new jump positions as start index for the random number ring
  for (i = 0; i < (nIterationAndLayer); i++)
  {
    idxRandomNumberStart[i] = floor(R::unif_rand() * (nIterationAndLayer));
  }
  PutRNGstate();

  uint32_t *sRange = R_Calloc(n * 2, uint32_t); 	// store for start and end positions of suitable ranges,
  // worst case amount of valid positions n * 2 (start & end for each filed, eq. checkerboard)
  size_t sRangeIdx; 								// index of "sRange", size_t
  uint32_t kgvLastRowIdx = ncol-1; 				    // ncol -1 left or right column as shortcut for ignoring surrounding cells
  uint32_t rangeStart;							    // index of the start index of a range in "sRange"
  bool isCounting;								    // flag beeing in suitable range or not

  // read origin d into xarrlT1
  // INFO: not required for xarrT0, xarrT0 will be initialized when copying xarrT1 to it
  // INFO: not required for xarrN3, xarrN3 is used in write before read (only written values will be used)
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for (i = 0; i < n; i++) {
    arrT1[i] = d[i]; // double -> int will always "floor"
  }

  // S layer loop (time periods)
  for (t = 0; t < nlay; t++) {

    // find suitable ranges
    sRangeIdx = 0;
    isCounting = false;
    for (i = ncol + 1; i < n - 1 - ncol; i++) { 					// ignore outer most cells of layer
      // start counting from first suitable postion
      if (!isCounting && (sPtr[i] > 0.0)) {
        // but dont start from first column OR the last
        if ((i % ncol == 0) || (i % ncol == kgvLastRowIdx)) {
          continue;
        }
        isCounting = true;
        sRange[sRangeIdx++] = i; 	// store start
        sRange[sRangeIdx] = i;		// store potential end
      }
      // end counting on last suitable position OR on second last column (to skip the last)
      else if (isCounting && (!(sPtr[i] > 0.0) || (i % ncol == kgvLastRowIdx))) {
        isCounting = false;
        sRange[sRangeIdx++] = i - 1; // store last suitable position as end of range
      }
    }
    // set the last suitable postion to the last element if it ends with only suitable positions
    if ((sRangeIdx % 2) != 0)
    {
      sRange[sRangeIdx++] = (n - 1) - (ncol + 1); // store last element - (last row - 1 col) as end of range
    }

    // iteration loop
    for (l = 0; l < it; l++) {

      // t0 <- t1
      memcpy(arrT0, arrT1, n * sizeof(uint8_t));

      // find and set occurrences within 3x3 neighborhood
      // parallel loop with thread private variables
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static) num_threads(nThreads) private(p, p1, p2, rangeStart, idxRandomRing) firstprivate(b)
      #endif
      for (i = 0; i < (sRangeIdx >> 1); i++) // >> division by 2, since only start positions should be used (eq. to increment "i" by 2)
      {
        rangeStart = i << 1; // get index of the start position in "sRange"
        // from start to end position
        for (p = sRange[rangeStart]; p <= sRange[rangeStart+1]; p++)
        {
          p1 = p+ncol;
          p2 = p-ncol;

          // get index of from buffer with random numbers,
          // wrap at the end of the buffer to use it as "ring buffer",
          // and prevent exceeding the number range of uint32_t by temporaily casting to uint64_t
          // further use of "idxRandomRing" requires the same modulo operations
          // index of random number must be unique through all layers, iterations, and positions in one layer
          // idxRandomRing = ((uint64_t)((p + (n * (t * it + l))) % randRingSize) * 11) % randRingSize;
          // idxRandomRing = ((((uint64_t)p + (n * (t * it + l)) + idxRandomNumberStart[currIterationAndLayer]) % randRingSize) * 11) % randRingSize;
          idxRandomRing = (
            (
              (uint64_t)idxRandomNumberStart[currIterationAndLayer] + // offset per iteration and layer
              ((uint64_t)p * 11) // 11 new random numbers per position (due to checks below, e.g. extinction, colonizaiton, corner events)
            ) % randRingSize // prevent index out of range
          );

          // extinction?
          if (arrT0[p]) {
            if (randRing[idxRandomRing] <= pe) { arrT1[p] = 0; }
          }

          // potential colonization only if focal cell uncolonized
          if (!arrT1[p]) {
            b = 0;    if (arrT0[p  ]) { b = (randRing[(idxRandomRing + 1) %randRingSize] <= sPtr[p]); }
            if (!b) { if (arrT0[p-1]) { b = (randRing[(idxRandomRing + 2) %randRingSize] <= sPtr[p]); } }
            if (!b) { if (arrT0[p+1]) { b = (randRing[(idxRandomRing + 3) %randRingSize] <= sPtr[p]); } }
            if (!b) { if (arrT0[p1 ]) { b = (randRing[(idxRandomRing + 4) %randRingSize] <= sPtr[p]); } }
            if (!b) { if (arrT0[p2 ]) { b = (randRing[(idxRandomRing + 5) %randRingSize] <= sPtr[p]); } }
          }
          // corner events
          if (!b) {
            if (randRing[(idxRandomRing + 6) %randRingSize] <= pc) {
              if (arrT0[p1-1]) { b = (randRing[(idxRandomRing + 7) %randRingSize] <= sPtr[p]); }
              if (!b) { if (arrT0[p1+1]) { b = (randRing[(idxRandomRing + 8) %randRingSize] <= sPtr[p]); } }
              if (!b) { if (arrT0[p2-1]) { b = (randRing[(idxRandomRing + 9) %randRingSize] <= sPtr[p]); } }
              if (!b) { if (arrT0[p2+1]) { b = (randRing[(idxRandomRing + 10) %randRingSize] <= sPtr[p]); } }
            }
          }
          //
          if ( b) { arrT1[p] = 1; }
        }
      }

      // calc (FOC=2; first colonization)
      if (ty == 2) {
        flocIdx = (t * it) + l + 1;
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nThreads)
        #endif
        for (i = 0; i < n; i++) { if ((arrT1[i]) && (ans[i] == 0)) { ans[i] = flocIdx; } }
      }

      // calc (LOC=3; last colonization)
      if (ty == 3) {
        flocIdx = (t*it)+l+1;
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nThreads)
        #endif
        for (i = 0; i < n; i++) { if (arrT1[i]) { ans[i] = flocIdx; } }
      }

      // calc (NOC=4; number of colonization events)
      if (ty == 4) {
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nThreads)
        #endif
        for (i = 0; i < n; i++) { if (arrT1[i]) { ans[i]++; } }
      }

      currIterationAndLayer++;
    }

    // set the first element of the new suitability layer as position 0, therefore sPtr[i] returns s[x*n+i] where n = layer size, x = current layer as index, i = the element of current layer
    sPtr += n;

  }

  // calc (DIS=1; final distribution)
  if (ty == 1) {
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for (i = 0; i < n; i++) { ans[i] = arrT1[i]; }
  }

  // signed?
  if (si == 1) {
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for (i = 0; i < n; i++) { if (arrT1[i] == 0) { ans[i] = -ans[i]; } }
  }

  // finally
  R_Free(arrT0);
  R_Free(arrT1);
  R_Free(occPos);
  R_Free(randRing);

  return ans;
}
