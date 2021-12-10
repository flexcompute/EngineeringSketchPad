
#include "Surreal/SurrealS.h"

#include "capsTypes.h"  // Bring in CAPS types
#include "miscUtils.h"

#include "vlmSpanSpace.h"

#include <vector>

#define PI        3.1415926535897931159979635
#define NINT(A)         (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))

template<class T>
void spacer( const T& N, const double pspace, T x[]) {
/*  modified from avl source sgutil.f
 *  computes only the spacing at the two ends of the segment
 *
 *...PURPOSE     TO CALCULATE A NORMALIZED (0<=X<=1) SPACING ARRAY.
 *
 *...INPUT       N      =  NUMBER OF DESIRED POINTS IN ARRAY.
 *               PSPACE =  SPACING PARAMETER (-3<=PSPACE<=3).
 *                         DEFINES POINT DISTRIBUTION
 *                         TO BE USED AS FOLLOWS:
 *                 PSPACE = 0  : EQUAL SPACING
 *                 PSPACE = 1  : COSINE SPACING.
 *                 PSPACE = 2  : SINE SPACING
 *                               (CONCENTRATING POINTS NEAR 0).
 *                 PSPACE = 3  : EQUAL SPACING.
 *
 *                 NEGATIVE VALUES OF PSPACE PRODUCE SPACING
 *                 WHICH IS REVERSED (AFFECTS ONLY SINE SPACING).
 *                 INTERMEDIATE VALUES OF PSPACE WILL PRODUCE
 *                 A SPACING WHICH IS A LINEAR COMBINATION
 *                 OF THE CORRESPONDING INTEGER VALUES.
 *
 *...OUTPUT      X      =  NORMALIZED SPACING ARRAY (0 <= X <= 1)
 *                         THE FIRST ELEMENT WILL ALWAYS BE  X(1) = 0.
 *                         THE LAST ELEMENT WILL ALWAYS BE   X(N) = 1.
 *
 */

  double pabs, pequ, pcos, psin;
  T frac, theta;
  int nabs;

  pabs = fabs(pspace);
  if (pabs > 3.0) pabs = 3;
  nabs = (int) pabs + 1;

  if (nabs == 1) {
    pequ   = 1.-pabs;
    pcos   = pabs;
    psin   = 0.;
  } else if (nabs == 2) {
    pequ   = 0.;
    pcos   = 2.-pabs;
    psin   = pabs-1.;
  } else {
    pequ   = pabs-2.;
    pcos   = 0.;
    psin   = 3.-pabs;
  }

  for (int i = 0; i < 2; i++) {
    if (i == 0) {
      frac = 1./(N-1.);
    } else {
      frac = (N-2.)/(N-1.);
    }
    theta =  frac * PI;
    if (pspace >= 0. ) x[i] =
                 pequ * frac
       +         pcos * ( 1. - cos ( theta )      ) / 2.
       +         psin * ( 1. - cos ( theta / 2. ) );
    if (pspace <= 0. ) x[i] =
                 pequ * frac
       +         pcos * ( 1. - cos ( theta )      ) / 2.
       +         psin * sin ( theta / 2. );
  }

  // change to a delta at the end
  x[1] = 1. - x[1];
}

// Compute auto spanwise panel spacing
extern "C"
int vlm_autoSpaceSpanPanels(int NspanTotal, int numSection, vlmSectionStruct vlmSection[])
{
    int    i, j, sectionIndex1, sectionIndex2;
    double distLE, distLETotal = 0, numSpanX;
    int numControl[2];
    int Nspan;
    int NspanMax, imax;
    int NspanMin, imin;

    int numSeg = numSection-1;
    
    // special case for just one segment (2 sections)
    if (numSeg == 1) {
        sectionIndex1 = vlmSection[0].sectionIndex;

        // use any specified counts
        if (vlmSection[sectionIndex1].Nspan >= 2)
          return CAPS_SUCCESS;

        // just set the total
        vlmSection[sectionIndex1].Nspan = NspanTotal;
        if (vlmSection[sectionIndex1].Sset == (int)false)
          vlmSection[sectionIndex1].Sspace = 1; // cosine

        return CAPS_SUCCESS;
    }

    // matrix/vectors for the newton solve
    std::vector<double> A(numSeg*numSeg, 0);
    std::vector<double> rhs(numSeg, 0);
    std::vector<double> x(numSeg, 0);
    std::vector<double> dx(numSeg, 0);

    std::vector<double> b(numSeg, 0); // lenght of each span section

    // go over all but the last section
    for (i = 0; i < numSection-1; i++) {

        // get the section indices
        sectionIndex1 = vlmSection[i  ].sectionIndex;
        sectionIndex2 = vlmSection[i+1].sectionIndex;

        // use the y-z distance between leading edge points to scale the number of spanwise points
        distLE = 0;
        for (j = 1; j < 3; j++) {
            distLE += pow(vlmSection[sectionIndex2].xyzLE[j] - vlmSection[sectionIndex1].xyzLE[j], 2);
        }
        distLE = sqrt(distLE);

        b[i] = distLE;
        distLETotal += distLE;

        // Modify the spacing parameter based on control surfaces present at each section
        if (vlmSection[sectionIndex1].Sset == (int)false) {

          numControl[0] = vlmSection[sectionIndex1].numControl;
          numControl[1] = vlmSection[sectionIndex2].numControl;

          // tips count as control surfaces as they need clustering as well
          // TODO: Deal with yduplicate
          if (i == 0           ) numControl[0] = 1;
          if (i == numSection-2) numControl[1] = 1;

               if (numControl[0] == 0 &&
                   numControl[1] == 0)
              vlmSection[sectionIndex1].Sspace = 0; // equal
          else if (numControl[0] > 0 &&
                   numControl[1] > 0)
              vlmSection[sectionIndex1].Sspace = 1; // cosine
          else if (numControl[0] > 0)
              vlmSection[sectionIndex1].Sspace = 2; // sine biast toward first section
          else if (numControl[1] > 0)
              vlmSection[sectionIndex1].Sspace = -2; // sine biast toward second section
        }
    }

    for (i = 0; i < numSeg; i++) {
        b[i] /= distLETotal;
        x[i] = b[i]*abs(NspanTotal);
    }

    // Using 5 newton iterations should be plenty without worrying about tolerances... (I hope)
    for (j = 0; j < 5; j++) {

        // reset the matrix
        for (i = 0; i < numSeg*numSeg; i++) A[i] = 0;
        for (i = 0; i < numSeg; i++) rhs[i] = 0;
        numSpanX = 0;
        for (i = 0; i < numSeg; i++) {
            if (x[i] < 3) x[i] = 3;
            numSpanX += x[i];
        }
        //printf("numSpanX = %f\n", numSpanX);

        // requite spacing on either side of a segment to be identical
        for (i = 0; i < numSeg; i++) {
            sectionIndex1 = vlmSection[i].sectionIndex;

            // use any specified counts
            if (vlmSection[sectionIndex1].Nspan >= 2) {
              A[numSeg*i + i] = 1;
              x[i] = vlmSection[sectionIndex1].Nspan;
              continue;
            }

            SurrealS<1> N = x[i]; N.deriv(0) = 1;

            SurrealS<1> dt[2];

            spacer( N, vlmSection[sectionIndex1].Sspace, dt);

            rhs[i]          = dt[1].value()*b[i];
            A[numSeg*i + i] = dt[1].deriv()*b[i];

            if (i > 0 && vlmSection[vlmSection[i-1].sectionIndex].Nspan < 2) {
                rhs[i-1]            -=  dt[0].value()*b[i];
                A[numSeg*(i-1) + i]  = -dt[0].deriv()*b[i];
            }
        }

        // The last equation requires the sum to add up to the total
        for (i = 0; i < numSeg; i++) {
            A[numSeg*(numSeg-1) + i] = 1;
        }
        rhs[numSeg-1] = numSpanX - abs(NspanTotal);

        // solve the linear system
        solveLU(numSeg, A.data(), rhs.data(), dx.data() );

        // update the solution
        for (i = 0; i < numSeg; i++) x[i] -= dx[i];
    }

    // set the number of spanwise points
    for (i = 0; i < numSection-1; i++) {

        // get the section indices
        sectionIndex1 = vlmSection[i].sectionIndex;

        vlmSection[sectionIndex1].Nspan = NINT(x[i]) > 2 ? NINT(x[i]) : 2;
    }


    // make sure the total adds up
    do {

        Nspan = 0;
        NspanMax = 0;
        NspanMin = NspanTotal;
        imax = 0;
        imin = 0;
        for (i = 0; i < numSection-1; i++) {
            // get the section indices
             sectionIndex1 = vlmSection[i].sectionIndex;

            if ( vlmSection[sectionIndex1].Nspan > NspanMax ) {
                NspanMax = vlmSection[sectionIndex1].Nspan;
                imax = sectionIndex1;
            }
            if ( vlmSection[sectionIndex1].Nspan < NspanMin ) {
                NspanMin = vlmSection[sectionIndex1].Nspan;
                imin = sectionIndex1;
            }

            Nspan += vlmSection[sectionIndex1].Nspan;
        }

        if (Nspan > NspanTotal) {
            vlmSection[imax].Nspan--;

            if (vlmSection[imax].Nspan == 1) {
                printf("Error: Insufficient spanwise sections! Increase numSpanTotal or numSpanPerSection!\n");
                return CAPS_BADVALUE;
            }
        }
        if (Nspan < NspanTotal) {
            vlmSection[imin].Nspan++;
        }

    } while (Nspan != NspanTotal);


    return CAPS_SUCCESS;
}
