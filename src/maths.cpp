/*==============================================================================*
 * MATHS
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/


#include "trace.h"
#include "maths.h"

std::mt19937 mt_gen(0);        /*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd() {return rnd_dist(mt_gen);}

/*computes determinant of a 4x4 matrix*/
double det4(double (*M)[4]) { TRACE_ME;
    double M0[3][3];
    double M1[3][3];
    double M2[3][3];
    double M3[3][3];

    for (int i=0;i<3;i++) {
        M0[i][0]=M[i+1][1];
        M0[i][1]=M[i+1][2];
        M0[i][2]=M[i+1][3];

        M1[i][0]=M[i+1][0];
        M1[i][1]=M[i+1][2];
        M1[i][2]=M[i+1][3];

        M2[i][0]=M[i+1][0];
        M2[i][1]=M[i+1][1];
        M2[i][2]=M[i+1][3];

        M3[i][0]=M[i+1][0];
        M3[i][1]=M[i+1][1];
        M3[i][2]=M[i+1][2];
    }

    return M[0][0]*det3(M0) -
           M[0][1]*det3(M1) +
           M[0][2]*det3(M2) -
           M[0][3]*det3(M3);
}

/*computes determinant of a 3x3 matrix*/
double det3(double (*M)[3]) { TRACE_ME;
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-
           M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+
           M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
}

/*helper functions for matrix math, y=A*x */
void matVecMultiply(double *y, double**A, double *x, int nu) { TRACE_ME;
    #pragma omp parallel for
    for (int i=0;i<nu;i++) {
        y[i] = 0;
        for (int j=0;j<nu;j++)
            y[i] += A[i][j]*x[j];
    }
}

/*computes y=v1-v2*/
void vecVecSubtract(double *y, double *v1, double *v2, int nu) {
    for (int i=0;i<nu;i++)
            y[i] = v1[i]-v2[i];
}


