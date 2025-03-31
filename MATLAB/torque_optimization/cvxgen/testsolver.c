/* Produced by CVXGEN, 2025-03-31 14:19:24 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
  return 0;
}
void load_default_data(void) {
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.P[0] = 1.5507979025745755;
  params.P[4] = 0;
  params.P[8] = 0;
  params.P[12] = 0;
  params.P[1] = 0;
  params.P[5] = 1.7081478226181048;
  params.P[9] = 0;
  params.P[13] = 0;
  params.P[2] = 0;
  params.P[6] = 0;
  params.P[10] = 1.2909047389129444;
  params.P[14] = 0;
  params.P[3] = 0;
  params.P[7] = 0;
  params.P[11] = 0;
  params.P[15] = 1.510827605197663;
  params.q[0] = 1.5717878173906188;
  params.q[1] = 1.5851723557337523;
  params.q[2] = -1.497658758144655;
  params.q[3] = -1.171028487447253;
  params.l[0] = -1.7941311867966805;
  params.l[1] = -0.23676062539745413;
  params.l[2] = -1.8804951564857322;
  params.l[3] = -0.17266710242115568;
  params.A[0] = 0.596576190459043;
  params.A[1] = -0.8860508694080989;
  params.A[2] = 0.7050196079205251;
  params.A[3] = 0.3634512696654033;
  params.A[4] = -1.9040724704913385;
  params.A[5] = 0.23541635196352795;
  params.A[6] = -0.9629902123701384;
  params.A[7] = -0.3395952119597214;
  params.A[8] = -0.865899672914725;
  params.A[9] = 0.7725516732519853;
  params.A[10] = -0.23818512931704205;
  params.A[11] = -1.372529046100147;
  params.A[12] = 0.17859607212737894;
  params.A[13] = 1.1212590580454682;
  params.A[14] = -0.774545870495281;
  params.A[15] = -1.1121684642712744;
  params.u[0] = -0.44811496977740495;
  params.u[1] = 1.7455345994417217;
  params.u[2] = 1.9039816898917352;
  params.u[3] = 0.6895347036512547;
}
