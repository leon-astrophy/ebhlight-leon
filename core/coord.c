/******************************************************************************
 *                                                                            *
 * COORD.C                                                                    *
 *                                                                            *
 * SIMULATION VOLUME COORDINATES                                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

/*      ASSUMING X3 SYMMETRY IN COORDINATES AND SPACETIME
 *      -- given the indices i,j and location in the cell, return with
 *         the values of X1,X2 there;
 *      -- the locations are defined by :
 *          -----------------------
 *          |                     |
 *          |                     |
 *          |FACE1   CENT         |
 *          |                     |
 *          |CORN    FACE2        |
 *          ----------------------
 *
 */

/**
 * @brief Fill the X matrix with the coordinates of a location on an active grid
 * 
 * @param i First zone index
 * @param j Second zone index
 * @param k Third zone index
 * @param loc Location within grid (0 = face1, 1 = face2, 2 = face3, 3 = corner, 4 = center)
 * @param X 1x4 matrix to fill with coords
 */
void coord(int i, int j, int k, int loc, double *X)
{
  X[0] = t;
	
  i += global_start[1];
  j += global_start[2];
  k += global_start[3];
  if (loc == FACE1) {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG)*dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else if (loc == FACE2) {
    X[1] = startx[1] + (i + 0.5 - NG)*dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else if (loc == FACE3) {
    X[1] = startx[1] + (i + 0.5 - NG)*dx[1];
    X[2] = startx[2] + (j + 0.5 - NG)*dx[2];
    X[3] = startx[3] + (k - NG) * dx[3];
  } else if (loc == CENT) {
    X[1] = startx[1] + (i + 0.5 - NG)*dx[1];
    X[2] = startx[2] + (j + 0.5 - NG)*dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else if (loc == CORN) {
    X[1] = startx[1] + (i - NG)*dx[1];
    X[2] = startx[2] + (j - NG)*dx[2];
    X[3] = startx[3] + (k - NG)*dx[3];
  } else {
    fprintf(stderr, "Invalid coordinate location!\n");
    exit(1);
  }
}

// Assumes Boyer-Lindquist coordinates
int bl_i_of_r(double r)
{
  double X1 = log(r);
  int ind = (int)((X1 - startx[1])/dx[1]);
  return ind;
}

double thG_of_X(const double X[NDIM])
{
  return M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
}

#if METRIC == MMKS
void thJ_of_X(const double X[NDIM], double *y, double* thJ)
{
  *y = 2*X[2] - 1.;
  *thJ = poly_norm*(*y)*(1. + pow((*y)/poly_xt,poly_alpha)/(poly_alpha+1.)) + 
    0.5*M_PI;
}
#endif

double r_of_X(const double X[NDIM])
{
  return exp(X[1]);
}

double th_of_X(const double X[NDIM])
{
  double thG = thG_of_X(X);
  
  #if METRIC == MMKS
  double y, thJ;
  thJ_of_X(X, &y, &thJ);
  return thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
  #else
  return thG;
  #endif
}

// Boyer-Lindquist coordinate of point X
void bl_coord(const double X[NDIM], double *r, double *th)
{
  *r = r_of_X(X);
  *th = th_of_X(X);

  // Avoid singularity at polar axis
  #if(COORDSINGFIX)
  if (fabs(*th) < SINGSMALL) {
    if ((*th) >= 0)
      *th = SINGSMALL;
    if ((*th) < 0)
      *th = -SINGSMALL;
  }
  if (fabs(M_PI - (*th)) < SINGSMALL) {
    if ((*th) >= M_PI)
      *th = M_PI + SINGSMALL;
    if ((*th) < M_PI)
      *th = M_PI - SINGSMALL;
  }
  #endif
}

// Cartesian coordinate Xcart = {t,x,y,z} of point X
/**
 * @brief Converts coordinates X to cartesian coordinates Xcart. If METRIC == CARTESIAN, Xcart == X
 * 
 * @param X Coordinates to convert (can be in cartesian, MKS or MMKS)
 * @param Xcart 1x4 matrix to fill with converted coordinates
 */
void cart_coord(const double X[NDIM], double Xcart[NDIM])
{
  Xcart[0] = X[0]; // t
  #if METRIC == MKS || METRIC == MMKS
  {
    double r,th,ph;
    bl_coord(X,&r,&th);
    ph = X[3];
    Xcart[1] = r*sin(th)*cos(ph); // r
    Xcart[2] = r*sin(th)*sin(ph); // th
    Xcart[3] = r*cos(th);         // ph
  }
  #else
  {
    Xcart[1] = X[1];
    Xcart[2] = X[2];
    Xcart[3] = X[3];
  }
  #endif // METRIC == MKS
}

// Jacobians.
// Jcon maps contravariant 4-vector from HARM coordinates
// to Boyer-Lindquist coordinates
// Jcov maps covariant 4-vector from HARM coordinates
// to Boyer-Lindquist coordinates
// CONVENTION: rows are top index, columns are bottom index
//             not required for tensor notation.
//             Is required for matrix notation.
void jac_harm_to_bl(const double X[NDIM],
  double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM])
{
  double r,th;
  bl_coord(X,&r,&th);

  double drdX1 = exp(X[1]);
  double dthGdX2 = M_PI + M_PI*(1 - hslope)*cos(2*M_PI*X[2]);
  double dthdX1,dthdX2;
  #if METRIC == MMKS
  {
    double thG = thG_of_X(X);
    double y, thJ;
    thJ_of_X(X,&y,&thJ);
    double dydX2 = 2.;
    double dthJdy = poly_norm*(1. + pow(y/poly_xt,poly_alpha));
    double dthJdX2 = dthJdy*dydX2;
    dthdX1 = -mks_smooth*(thG-thJ)*exp(mks_smooth*(startx[1]-X[1]));
    dthdX2 = dthGdX2 + exp(mks_smooth*(startx[1]-X[1]))*(dthGdX2 - dthJdX2);
  }
  #else
  {
    dthdX1 = 0.0;
    dthdX2 = dthGdX2;
  }
  #endif

  // First set all values to zero. Then we'll set nonzero values
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Jcon[mu][nu] = 0.;
      Jcov[mu][nu] = 0.;
    }
  }
  // Jcon
  Jcon[0][0] = 1.;     // t
  Jcon[1][1] = drdX1;  // r
  Jcon[2][1] = dthdX1; // th
  Jcon[2][2] = dthdX2; 
  Jcon[3][3] = 1.;     // phi
  // Jcov
  Jcov[0][0] = 1.;                             // t
  Jcov[1][1] = 1./(drdX1 + SMALL);             // r
  Jcov[2][1] = -dthdX1/(drdX1*dthdX2 + SMALL); // th
  Jcov[2][2] = 1./(dthdX2 + SMALL);
  Jcov[3][3] = 1.;
}

// Jcon maps contravariant 4-vector from BL coordinates to Minkowski
// Jcov maps a covariant one from BL to Minkowski
void jac_bl_to_cart(const double X[NDIM],
  double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM])
{
  double r,th,ph;
  double Xcart[NDIM];

  bl_coord(X,&r,&th);
  ph = X[3];
  cart_coord(X,Xcart);

  // First set all values to zero. Then we'll set nonzero values
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Jcon[mu][nu] = 0.;
      Jcov[mu][nu] = 0.;
    }
  }
  // Jcon
  Jcon[0][0] = 1.;                 // t
  Jcon[1][1] = sin(th)*cos(ph);    // x
  Jcon[1][2] = r*cos(th)*cos(ph);
  Jcon[1][3] = -r*sin(th)*sin(ph);
  Jcon[2][1] = sin(th)*sin(ph);    // y
  Jcon[2][2] = r*cos(th)*sin(ph);
  Jcon[2][3] = r*sin(th)*cos(ph);
  Jcon[3][1] = cos(th);            // z
  Jcon[3][2] = -r*sin(th);
  // Jcov
  Jcov[0][0] = 1.;
  Jcov[1][1] = cos(ph)*sin(th);
  Jcov[1][2] = sin(ph)*sin(th);
  Jcov[1][3] = cos(th);
  Jcov[2][1] = cos(ph)*cos(th)/(r+SMALL);
  Jcov[2][2] = cos(th)*sin(ph)/(r+SMALL);
  Jcov[2][3] = -sin(th)/(r+SMALL);
  Jcov[3][1] = -sin(ph)/(r*sin(th)+SMALL);
  Jcov[3][2] = cos(ph)/(r*sin(th)+SMALL);
}

void jac_harm_to_cart(const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM])
{
  #if METRIC == MKS || METRIC == MMKS
  {
    double J_h2bl_cov[NDIM][NDIM], J_h2bl_con[NDIM][NDIM];
    double J_bl2c_cov[NDIM][NDIM], J_bl2c_con[NDIM][NDIM];
    jac_harm_to_bl(X, J_h2bl_cov, J_h2bl_con);
    jac_bl_to_cart(X, J_bl2c_cov, J_bl2c_con);
    
    for (int mupp = 0; mupp < NDIM; mupp++) {
      for (int mu = 0; mu < NDIM; mu++) {
        Jcon[mupp][mu] = 0.0;
        Jcov[mu][mupp] = 0.0;
        for (int mup = 0; mup < NDIM; mup++) {
          Jcon[mupp][mu] += J_bl2c_con[mupp][mup]*J_h2bl_con[mup][mu];
          Jcov[mu][mupp] += J_h2bl_cov[mu][mup]*J_bl2c_cov[mup][mupp];
        }
      }
    }
  }
  #else
  {
    DLOOP2 Jcon[mu][nu] = Jcov[mu][nu] = (mu == nu ? 1 : 0);
  }
  #endif // METRIC
}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  memset(dxdX, 0, NDIM*NDIM*sizeof(double));
  #if METRIC == CARTESIAN
  for (int mu = 0; mu < NDIM; mu++) {
    dxdX[mu][mu] = 1.;
  }
  #elif METRIC == MKS
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
  dxdX[3][3] = 1.;
  #elif METRIC == MMKS
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
    M_PI/2. - 
    M_PI*X[2] + 
    poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) - 
    1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
    );
  dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) + 
    exp(mks_smooth*(startx[1]-X[1]))*(
      -M_PI + 
      2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) + 
      (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) - 
      (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
      );
  dxdX[3][3] = 1.;
  #endif
}

/**
 * @brief Set the covariant g object at given coords
 * 
 * @param X 1x4 matrix of coords
 * @param gcov 4x4 matrix to set to gcov
 */
void set_gcov(double X[NDIM], double gcov[NDIM][NDIM])
{
  memset(gcov, 0, NDIM*NDIM*sizeof(double));

  #if METRIC == CARTESIAN
  gcov[0][0] = -1.;
  for (int j = 1; j < NDIM; j++) {
    gcov[j][j] = 1.;
  }
  #elif METRIC == MKS || METRIC == MMKS
  double sth, cth, s2, rho2;
  double r, th;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  double gcov_ks[NDIM][NDIM];
  memcpy(gcov_ks, gcov, NDIM*NDIM*sizeof(double));
  memset(gcov, 0, NDIM*NDIM*sizeof(double));

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int lam = 0; lam < NDIM; lam++) {
        for (int kap = 0; kap < NDIM; kap++) {
          gcov[mu][nu] += gcov_ks[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
        }
      }
    }
  }
  #endif // METRIC
}

/**
 * @brief Establish starting, stopping coords and step size (startx, startx_rad, startx_proc, stopx, stopx_rad, stopx_proc, dx)
 * 
 */
void set_points()
{
  #if METRIC == CARTESIAN
  startx[1] = x1Min;
  startx[2] = x2Min;
  startx[3] = x3Min;
  dx[1] = (x1Max - x1Min)/N1TOT;
  dx[2] = (x2Max - x2Min)/N2TOT;
  dx[3] = (x3Max - x3Min)/N3TOT;
  #if RADIATION
  startx_rad[1] = startx[1];
  startx_rad[2] = startx[2];
  startx_rad[3] = startx[3];
  stopx_rad[1] = startx_rad[1] + N1TOT*dx[1];
  stopx_rad[2] = startx_rad[2] + N2TOT*dx[2];
  stopx_rad[3] = startx_rad[3] + N3TOT*dx[3];
  #endif
  #elif METRIC == MKS || METRIC == MMKS 
  // Set Rin such that we have 5 zones completely inside the event horizon
  Reh = 1. + sqrt(1. - a*a);
  double z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 + a, 1./3.) + pow(1 - a, 1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  Risco = 3 + z2 - sqrt((3 - z1)*(3 + z1 + 2*z2));
  Rin = exp((N1TOT*log(Reh)/5.5 - log(Rout))/(1. + N1TOT/5.5));

  // Inner radius outside event horizon specifically for thermalizing atmosphere
  if (!strcmp(PROBLEM_NAME, "bhtherm")) {
    Rin = 3.5;
    #if RADIATION
    Rout_rad = 20.;
    #endif
  }

  startx[1] = log(Rin);
  startx[2] = 0.;
  startx[3] = 0.;
  dx[1] = log(Rout/Rin)/N1TOT;
  dx[2] = 1./N2TOT;
  dx[3] = 2.*M_PI/N3TOT;

  #if RADIATION
    startx_rad[1] = log(Reh);
    if (!strcmp(PROBLEM_NAME, "bhtherm")) {
      startx_rad[1] = log(Rin);
      startx[2] = 0.05;
      dx[2] = 0.9/N2TOT;
    }
    startx_rad[2] = startx[2];
    startx_rad[3] = startx[3];
    stopx_rad[1] = log(Rout_rad);
    stopx_rad[2] = startx[2] + N2TOT*dx[2];
    stopx_rad[3] = startx[3] + N3TOT*dx[3];
  #endif
  #if METRIC == MMKS
    poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*
                           1./pow(poly_xt, poly_alpha));
  #endif 
  #endif // MKS

  startx_proc[1] = startx[1] + global_start[1]*dx[1];
  startx_proc[2] = startx[2] + global_start[2]*dx[2];
  startx_proc[3] = startx[3] + global_start[3]*dx[3];
  stopx_proc[1] = startx_proc[1] + dx[1]*N1;
  stopx_proc[2] = startx_proc[2] + dx[2]*N2;
  stopx_proc[3] = startx_proc[3] + dx[3]*N3;

  stopx[1] = startx[1] + dx[1]*N1TOT;
  stopx[2] = startx[2] + dx[2]*N2TOT;
  stopx[3] = startx[3] + dx[3]*N3TOT;
}

/**
 * @brief Set boundary and steps, then loop through active zones setting g matrix, g, alpha and connection coefficient
 * 
 */
void set_grid()
{
  double X[NDIM];

  // Set up boundaries, steps in coordinate grid
  set_points();
  dV = dx[1]*dx[2]*dx[3];

  for (int mu = 0; mu < NDIM; mu++) X[mu] = 0.;

  #if RADIATION
  dt_light_min = 1./SMALL;
  #endif 

  ISLOOP(-NG, N1 - 1 + NG) {
    JSLOOP(-NG, N2 - 1 + NG) {
      LOCLOOP {
        coord(i, j, 0, loc, X);
        set_gcov(X, ggeom[i][j][loc].gcov);
        ggeom[i][j][loc].g = gcon_func(ggeom[i][j][loc].gcov, ggeom[i][j][loc].gcon);
        ggeom[i][j][loc].alpha = 1.0 / sqrt(-ggeom[i][j][loc].gcon[0][0]);
      }
       
      // Only required in zone center
      conn_func(X, &ggeom[i][j][CENT], conn[i][j]);
      
      #if RADIATION
      // Set minimum light crossing time for each zone
      dt_light[i][j] = 1.e30;
      double light_phase_speed = SMALL;
      double dt_light_local = 0.; 

      for (int mu = 1; mu < NDIM; mu++) {
        if(pow(ggeom[i][j][CENT].gcon[0][mu], 2.) -
           ggeom[i][j][CENT].gcon[mu][mu]*ggeom[i][j][CENT].gcon[0][0] >= 0.) {
          double cplus = fabs((-ggeom[i][j][CENT].gcon[0][mu] +
            sqrt(pow(ggeom[i][j][CENT].gcon[0][mu], 2.) -
            ggeom[i][j][CENT].gcon[mu][mu]*ggeom[i][j][CENT].gcon[0][0]))/
            (ggeom[i][j][CENT].gcon[0][0]));
          double cminus = fabs((-ggeom[i][j][CENT].gcon[0][mu] -
            sqrt(pow(ggeom[i][j][CENT].gcon[0][mu], 2.) -
            ggeom[i][j][CENT].gcon[mu][mu]*ggeom[i][j][CENT].gcon[0][0]))/
            (ggeom[i][j][CENT].gcon[0][0])) ;
          light_phase_speed= MY_MAX(cplus,cminus);
        } else {
          light_phase_speed = SMALL;
        }

        dt_light_local += 1./(dx[mu]/light_phase_speed);

        if (dx[mu]/light_phase_speed < dt_light[i][j]) {
          dt_light[i][j] = dx[mu]/light_phase_speed;
        }
      }
      dt_light_local = 1./dt_light_local;
      if (dt_light_local < dt_light_min) dt_light_min = dt_light_local;
      #endif // RADIATION
    } // JSLOOP
  } // ISLOOP

  #if RADIATION
  dt_light_min = mpi_min(dt_light_min);
  #endif
}

/**
 * @brief Loop through all active zones, setting each primitive and flux to 0. Also 0s the pflag and fail_save flags in each zone
 * 
 */
void zero_arrays()
{
  ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG, -NG, N3-1 + NG) {
    PLOOP {
      P[i][j][k][ip] = 0.;
      Ph[i][j][k][ip] = 0.;
      F1[i][j][k][ip] = 0.;
      F2[i][j][k][ip] = 0.;
      F3[i][j][k][ip] = 0.;
    }
    pflag[i][j][k] = 0;
    fail_save[i][j][k] = 0;
  } // ZSLOOP
}

