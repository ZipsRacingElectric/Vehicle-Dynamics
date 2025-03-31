/* Produced by CVXGEN, 2025-03-31 14:19:24 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-params.A[0])-rhs[1]*(-params.A[4])-rhs[2]*(-params.A[8])-rhs[3]*(-params.A[12]);
  lhs[1] = -rhs[0]*(-params.A[1])-rhs[1]*(-params.A[5])-rhs[2]*(-params.A[9])-rhs[3]*(-params.A[13]);
  lhs[2] = -rhs[0]*(-params.A[2])-rhs[1]*(-params.A[6])-rhs[2]*(-params.A[10])-rhs[3]*(-params.A[14]);
  lhs[3] = -rhs[0]*(-params.A[3])-rhs[1]*(-params.A[7])-rhs[2]*(-params.A[11])-rhs[3]*(-params.A[15]);
  lhs[4] = -rhs[0]*(params.A[0])-rhs[1]*(params.A[4])-rhs[2]*(params.A[8])-rhs[3]*(params.A[12]);
  lhs[5] = -rhs[0]*(params.A[1])-rhs[1]*(params.A[5])-rhs[2]*(params.A[9])-rhs[3]*(params.A[13]);
  lhs[6] = -rhs[0]*(params.A[2])-rhs[1]*(params.A[6])-rhs[2]*(params.A[10])-rhs[3]*(params.A[14]);
  lhs[7] = -rhs[0]*(params.A[3])-rhs[1]*(params.A[7])-rhs[2]*(params.A[11])-rhs[3]*(params.A[15]);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-params.A[0])-rhs[1]*(-params.A[1])-rhs[2]*(-params.A[2])-rhs[3]*(-params.A[3])-rhs[4]*(params.A[0])-rhs[5]*(params.A[1])-rhs[6]*(params.A[2])-rhs[7]*(params.A[3]);
  lhs[1] = -rhs[0]*(-params.A[4])-rhs[1]*(-params.A[5])-rhs[2]*(-params.A[6])-rhs[3]*(-params.A[7])-rhs[4]*(params.A[4])-rhs[5]*(params.A[5])-rhs[6]*(params.A[6])-rhs[7]*(params.A[7]);
  lhs[2] = -rhs[0]*(-params.A[8])-rhs[1]*(-params.A[9])-rhs[2]*(-params.A[10])-rhs[3]*(-params.A[11])-rhs[4]*(params.A[8])-rhs[5]*(params.A[9])-rhs[6]*(params.A[10])-rhs[7]*(params.A[11]);
  lhs[3] = -rhs[0]*(-params.A[12])-rhs[1]*(-params.A[13])-rhs[2]*(-params.A[14])-rhs[3]*(-params.A[15])-rhs[4]*(params.A[12])-rhs[5]*(params.A[13])-rhs[6]*(params.A[14])-rhs[7]*(params.A[15]);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(params.P[0])+rhs[1]*(params.P[4])+rhs[2]*(params.P[8])+rhs[3]*(params.P[12]);
  lhs[1] = rhs[0]*(params.P[1])+rhs[1]*(params.P[5])+rhs[2]*(params.P[9])+rhs[3]*(params.P[13]);
  lhs[2] = rhs[0]*(params.P[2])+rhs[1]*(params.P[6])+rhs[2]*(params.P[10])+rhs[3]*(params.P[14]);
  lhs[3] = rhs[0]*(params.P[3])+rhs[1]*(params.P[7])+rhs[2]*(params.P[11])+rhs[3]*(params.P[15]);
}
void fillq(void) {
  work.q[0] = params.q[0];
  work.q[1] = params.q[1];
  work.q[2] = params.q[2];
  work.q[3] = params.q[3];
}
void fillh(void) {
  work.h[0] = -params.l[0];
  work.h[1] = -params.l[1];
  work.h[2] = -params.l[2];
  work.h[3] = -params.l[3];
  work.h[4] = params.u[0];
  work.h[5] = params.u[1];
  work.h[6] = params.u[2];
  work.h[7] = params.u[3];
}
void fillb(void) {
}
void pre_ops(void) {
}
