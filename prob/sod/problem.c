/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOD SHOCKTUBE                                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "io.h"

static double plaw;
static double tscale;
void set_problem_params()
{
  set_param("plaw", &plaw);
  set_param("tscale", &tscale);
}
void save_problem_params()
{
  WRITE_HDR(plaw, TYPE_DBL);
  WRITE_HDR(tscale, TYPE_DBL);
}

void init_prob()
{
	double X[NDIM];

  // Make problem nonrelativistic
  double tscale = 1.e-2;
  tf /= tscale;
  dt /= tscale;
  DTd /= tscale;
  DTl /= tscale;

  ZLOOP {
    coord(i, j, k, CENT, X);

    PLOOP P[i][j][k][ip] = 0.;

    P[i][j][k][RHO] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;

    double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
    P[i][j][k][UU] = pgas/(gam - 1.);

    #if NONTHERMAL
    NTELOOP P[i][j][k][ip] = (X[1] < 0.5 || X[1] > 1.5) ? 1e-10 : 0.125e-10;
    #endif
  } // ZLOOP

  // Rescale to make problem nonrelativistic
  ZLOOP {
    P[i][j][k][UU] *= tscale*tscale;
    P[i][j][k][U1] *= tscale;
    P[i][j][k][U2] *= tscale;
    P[i][j][k][U3] *= tscale;
    P[i][j][k][B1] *= tscale;
    P[i][j][k][B2] *= tscale;
    P[i][j][k][B3] *= tscale;
  }
}

