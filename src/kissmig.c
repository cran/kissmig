/* Michael Nobis, May 2014 */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>    // R RNG
#include "Rdefines.h"
#include "R_ext/Rdynload.h"


SEXP kissmig_c(SEXP d, SEXP s, SEXP dim, SEXP it, SEXP pe, SEXP pc, SEXP ty, SEXP si) {

	SEXP ans;
	SEXP arrT0;
	SEXP arrT1;
	SEXP arrN3;

	int i, j, t, l, p, ps, p1, p2, nrow, ncol, nlay, n, b;
	int *xarrT0, *xarrT1, *xarrN3;
	double *xans, *xd, *xs;

	int xit = INTEGER(it)[0];
	int xty = INTEGER(ty)[0];
	int xsi = INTEGER(si)[0];
	double xpe = REAL(pe)[0];
	double xpc = REAL(pc)[0];

	nrow = INTEGER(dim)[0];
	ncol = INTEGER(dim)[1];
	nlay = INTEGER(dim)[2];
	n = nrow * ncol;

	PROTECT( d     = coerceVector(d, REALSXP));
	PROTECT( s     = coerceVector(s, REALSXP));
	PROTECT( ans   = allocVector(REALSXP, n) );
	PROTECT( arrT0 = allocVector(LGLSXP , n) );
	PROTECT( arrT1 = allocVector(LGLSXP , n) );
	PROTECT( arrN3 = allocVector(LGLSXP , n) );

	xd     = REAL(d);
	xs     = REAL(s);
	xans   = REAL(ans);
	xarrT0 = INTEGER(arrT0);
	xarrT1 = INTEGER(arrT1);
	xarrN3 = INTEGER(arrN3);

  	// R RNG
 	GetRNGstate();

	// read origin d into xarrlT1 and initialize xarrT0, xarrN3x3
    for (i = 0; i < n; i++) {
		xarrT1[i] = (xd[i] == 1);
		xarrT0[i] = 0;
		xarrN3[i] = 0;
		xans[i]   = 0;
	}

	// S layer loop (time periods)
	for (t = 0; t < nlay; t++) {

		// iteration loop
		for (l = 0; l < xit; l++) {

			// t0 <- t1
			for (i = 0; i < n; i++) { xarrT0[i] = xarrT1[i]; }

			// occurrences within 3x3 neighborhood (opt?)
			p = ncol-1;
			for (i = 1; i < nrow-1; i++) {
				p++;
				for (j = 1; j < ncol-1; j++) {
					p++;
					p1 = p+ncol;
					p2 = p-ncol;
					xarrN3[p] =
						(xarrT0[p ]) || (xarrT0[p -1]) || (xarrT0[p +1]) ||
						(xarrT0[p1]) || (xarrT0[p1-1]) || (xarrT0[p1+1]) ||
						(xarrT0[p2]) || (xarrT0[p2-1]) || (xarrT0[p2+1]);
				}
				p++;
			}

			// new distribution at T1
			p  = ncol-1;
			ps = (t*n)+p;
			for (i = 1; i < nrow-1; i++) {
				p++; ps++;
				for (j = 1; j < ncol-1; j++) {
					p++; ps++;
					p1 = p+ncol;
					p2 = p-ncol;

					// potential colonization only if occurrences within 3x3 neighborhood
					if (xarrN3[p]) {

						// extinction?
						if (xarrT1[p]) { if (unif_rand()<=xpe) {xarrT1[p] = 0; } }

						// potential colonization only if focal cell uncolonized
						b = 0;
						if (!xarrT1[p]) {
							          if (xarrT0[p  ]) { b = (unif_rand() <= xs[ps]); }
							if (!b) { if (xarrT0[p-1]) { b = (unif_rand() <= xs[ps]); } }
							if (!b) { if (xarrT0[p+1]) { b = (unif_rand() <= xs[ps]); } }
							if (!b) { if (xarrT0[p1 ]) { b = (unif_rand() <= xs[ps]); } }
							if (!b) { if (xarrT0[p2 ]) { b = (unif_rand() <= xs[ps]); } }
						}
						// corner events
						if (!b) {
							if (unif_rand() <= xpc) {
										  if (xarrT0[p1-1]) { b = (unif_rand() <= xs[ps]); }
								if (!b) { if (xarrT0[p1+1]) { b = (unif_rand() <= xs[ps]); } }
								if (!b) { if (xarrT0[p2-1]) { b = (unif_rand() <= xs[ps]); } }
								if (!b) { if (xarrT0[p2+1]) { b = (unif_rand() <= xs[ps]); } }
							}
						}
						//
						if ( b) { xarrT1[p] = 1; }

					}
				}
				p++; ps++;

			}

			// calc (FCO=2; first colonization)
			if (xty == 2) { for (i = 0; i < n; i++) { if ((xarrT1[i]) && (xans[i] == 0)) { xans[i] = (t*xit)+l+1; } } }

			// calc (LCO=3; last colonization)
			if (xty == 3) { for (i = 0; i < n; i++) { if (xarrT1[i]) { xans[i] = (t*xit)+l+1; } } }

			// calc (NCO=4; number of colonization events)
			if (xty == 4) { for (i = 0; i < n; i++) { if (xarrT1[i]) { xans[i]++; } } }

		}

	}

	// calc (DIS=1; final distribution)
	if (xty == 1) { for (i = 0; i < n; i++) { xans[i] = xarrT1[i]; } }

	// signed?
	if (xsi == 1) { for (i = 0; i < n; i++) { if (xarrT1[i] == 0) { xans[i] = -xans[i]; } } }

	// finally
 	PutRNGstate();
	UNPROTECT(6);

	return(ans);
}

