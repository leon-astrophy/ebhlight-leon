/******************************************************************************
 *                                                                            *
 * STEP.C                                                                     *
 *                                                                            *
 * ADVANCES FLUID QUANTITIES BY ONE TIMESTEP                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt,
  grid_prim_type Pf, int stage);
double fluxcalc(grid_prim_type Pr);
void flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
void lr_to_flux(double p_l[NVAR], double p_r[NVAR], struct of_geom *geom,
  int dir, double Flux[NVAR], double *maxSpeed);
#if RADIATION
void apply_rad_force(grid_prim_type Pr, double Dt);
#endif

void step()
{
  double ndt;
  dtsave = dt;

  // Need both P_n and P_n+1 to calculate current
  ZSLOOP(-NG,N1-1+NG,-NG,N2-1+NG,-NG,N3-1+NG) {
    #ifdef CONST_INJECTION
      inject_nonthermal(P[i][j][k], CONST_INJECTION, dt);
    #endif
    
    PLOOP {
      Psave[i][j][k][ip] = P[i][j][k][ip];
    }
  }

  // Predictor step - advance half a timestep
  ndt = advance(P, P, 0.5*dt, Ph, 0);
  #if NONTHERMAL && !defined(SKIP_ADIAB)
  nonthermal_adiab(P, Ph, 0.5*dt);
  #endif
  #if ELECTRONS
  heat_electrons(P, P, Ph, 0.5*dt);
  #endif
  #if NONTHERMAL && !defined(SKIP_COOLING)
  cool_nonthermal(P, Ph, 0.5*dt);
  #endif
  fixup(Ph);
  fixup_utoprim(Ph);
  #if ELECTRONS
  fixup_electrons(Ph);
  #if RADIATION
  coulomb(P, P, Ph, 0.5*dt);
  fixup_electrons(Ph);
  #endif
  #endif
  bound_prim(Ph);

  if (mpi_io_proc())
  // Radiation step
  #if RADIATION
  make_superphotons(Ph, t, dt);
  if (mpi_io_proc())
    //fprintf(stdout, "Made superphotons\n");
  push_superphotons(dt);
  if (mpi_io_proc())
    //fprintf(stdout, "Pushed superphotons\n");
  interact(Ph, t, dt);
  if (mpi_io_proc())
    //fprintf(stdout, "Interacted\n");
  bound_superphotons(t, dt);
  if (mpi_io_proc())
    //fprintf(stdout, "Bound\n");
  #endif
  // Corrector step
  ndt = advance(P, Ph, dt, P, 1);
  #if NONTHERMAL && !defined(SKIP_ADIAB)
  nonthermal_adiab(Psave, P, dt);
  #endif
  #if ELECTRONS
  heat_electrons(P, Ph, P, dt);
  #endif
  #if NONTHERMAL && !defined(SKIP_COOLING)
  cool_nonthermal(P, P, dt);
  #endif
  fixup(P);
  fixup_utoprim(P);
  #if ELECTRONS
  fixup_electrons(P);
  #if RADIATION
  coulomb(P, Ph, P, dt);
  fixup_electrons(P);
  #endif
  #endif
  bound_prim(P);

  #if RADIATION

  // Apply radiation four-force to fluid
  apply_rad_force(P, dt);
  fixup(P);
  fixup_utoprim(P);
  #if ELECTRONS
  #if NONTHERMAL
  apply_rad_force_e(Ph, P, radG_e, dt);
  #else
  apply_rad_force_e(Ph, P, radG, dt);
  #endif
  fixup_electrons(P);
  #endif
  bound_prim(P);

  memset((void*)&radG[0][0][0][0], 0,
    (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*NDIM*sizeof(double));

  #if NONTHERMAL
  memset((void*)&radG_e[0][0][0][0], 0,
    (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*NDIM*sizeof(double));
  #endif
  #endif // RADIATION

  // Increment time
  t += dt;
  
  #if RADIATION
  ndt = cour*dt_light_min; 
  #endif

  // Set next timestep
  if (ndt > SAFE * dt) {
    ndt = SAFE * dt;
  }

  dt = ndt;

  // Don't step beyond end of run
  if (t + dt > tf) {                                                         
    dt = tf - t;                                                             
  }
  
  dt = mpi_min(dt);
}

/**
 * @brief Advances the Primitive matrix forward by a timestep dt. Does F(Pb) -> U(F) -> Pf(U)
 * 
 * @param Pi Current integer step primitives
 * @param Pb Previous half-step primitives
 * @param Dt Timestep over which to evolve Pi to Pf
 * @param Pf Matrix to store final primitive values
 * @param stage Unused
 * @return double ndt = 1/(1/ndt1 + 1/ndt2 + 1/ndt3)
 */
double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt,
  grid_prim_type Pf, int stage)
{
  double ndt, U[NVAR], dU[NVAR];
  struct of_state qi;

  #pragma omp parallel for collapse(3)
  ZLOOP PLOOP Pf[i][j][k][ip] = Pi[i][j][k][ip];

  timer_start(TIMER_FLUXCALC);
  ndt = fluxcalc(Pb);

  #if METRIC == MKS || METRIC == MMKS
  fix_flux(F1, F2, F3);
  #endif

  flux_ct(F1, F2, F3);
  timer_stop(TIMER_FLUXCALC);

  // Evaluate diagnostics based on fluxes
  timer_start(TIMER_DIAG);
  diag_flux(F1, F2, F3);
  timer_stop(TIMER_DIAG);

  // Update Pi to Pf
  timer_start(TIMER_UPDATE);
  #pragma omp parallel for schedule(guided) private (dU, qi, U ) collapse(3)
  ZLOOP {
    source(Pb[i][j][k], &(ggeom[i][j][CENT]), i, j, dU, Dt);
    get_state(Pi[i][j][k], &(ggeom[i][j][CENT]), &qi);
    primtoflux(Pi[i][j][k], &qi, 0, &(ggeom[i][j][CENT]), U);

    PLOOP {
      U[ip] += Dt*(
       - (F1[i + 1][j][k][ip] - F1[i][j][k][ip]) / dx[1]
       - (F2[i][j + 1][k][ip] - F2[i][j][k][ip]) / dx[2]
       - (F3[i][j][k + 1][ip] - F3[i][j][k][ip]) / dx[3]
       + dU[ip]);
    }

    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pf[i][j][k]);
    if(pflag[i][j][k])
      fail_save[i][j][k] = 1;
  } // ZLOOP
  timer_stop(TIMER_UPDATE);

  return (ndt);
}

#if RADIATION
void apply_rad_force(grid_prim_type Pr, double Dt)
{
  double U[NVAR];
  struct of_state q;
  
  timer_start(TIMER_UPDATE);

  #pragma omp parallel for schedule(guided) private (q, U) collapse(3)
  ZLOOP {
    // Store primitive variables before cooling for supercooling diagnostics
    PLOOP psupersave[i][j][k][ip] = Pr[i][j][k][ip];

    get_state(Pr[i][j][k], &(ggeom[i][j][CENT]), &q);
    primtoflux(Pr[i][j][k], &q, 0, &(ggeom[i][j][CENT]), U);

    for (int ip = 1; ip < 5; ip++) {
      U[ip] += Dt*radG[i][j][k][ip-1];
      radG_prev[i][j][k][ip-1] = radG[i][j][k][ip-1];
    }

    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pr[i][j][k]);
    
    if(pflag[i][j][k]) {
      fail_save[i][j][k] = 1;
    }
  } // ZLOOP
  timer_stop(TIMER_UPDATE);
}
#endif

/**
 * @brief Find the fluxes from the primitives and populate F1, F2, and F3 with flux in each spatial dimension
 * 
 * @param Pr Primitives at each gridpoint
 * @return double ndt = 1/(1/ndt1 + 1/ndt2 + 1/ndt3)
 */
double fluxcalc(grid_prim_type Pr)
{
  double P_l[NMAX+2*NG][NVAR], P_r[NMAX+2*NG][NVAR];
  double cij, cmax1, cmax2, cmax3;
  double Ptmp[NMAX+2*NG][NVAR];

  cmax1 = cmax2 = cmax3 = 0.;

  #pragma omp parallel private(Ptmp, P_l, P_r, cij) \
    reduction(max:cmax1) reduction(max:cmax2) reduction(max:cmax3)
  {
    #pragma omp for collapse(2) nowait
    JSLOOP(-1,N2) {
      KSLOOP(-1,N3) {
        ISLOOP(-NG,N1-1+NG) PLOOP Ptmp[i][ip] = Pr[i][j][k][ip];
        
        reconstruct(Ptmp, N1, P_l, P_r);

        ISLOOP(0,N1) {
          lr_to_flux(P_r[i-1], P_l[i], &(ggeom[i][j][FACE1]), 1, F1[i][j][k], &cij);
          cmax1 = (cij > cmax1 ? cij : cmax1);
        } // ISLOOP
      } // KSLOOP
    } // JSLOOP

    #pragma omp for collapse(2) nowait
    ISLOOP(-1,N1) {
      KSLOOP(-1,N3) {

        JSLOOP(-NG,N2-1+NG) PLOOP Ptmp[j][ip] = Pr[i][j][k][ip];

        reconstruct(Ptmp, N2, P_l, P_r);

        JSLOOP(0,N2) {
          lr_to_flux(P_r[j-1], P_l[j], &(ggeom[i][j][FACE2]), 2, F2[i][j][k], &cij);
          cmax2 = (cij > cmax2 ? cij : cmax2);
        } // JSLOOP
      } // KSLOOP
    } // ISLOOP

    #pragma omp for collapse(2)
    ISLOOP(-1,N1) {
      JSLOOP(-1,N2) {

        KSLOOP(-NG,N3-1+NG) PLOOP Ptmp[k][ip] = Pr[i][j][k][ip];

        reconstruct(Ptmp, N3, P_l, P_r);

        KSLOOP(0,N3) {
          lr_to_flux(P_r[k-1], P_l[k], &(ggeom[i][j][FACE3]), 3, F3[i][j][k], &cij);
          cmax3 = (cij > cmax3 ? cij : cmax3);
        } // KSLOOP
      } // JSLOOP
    } // ISLOOP
  } // omp parallel

  // Otherwise timestep changes with MPI!
  cmax1 = mpi_max(cmax1);
  cmax2 = mpi_max(cmax2);
  cmax3 = mpi_max(cmax3);

  double ndt1 = cour * dx[1] / cmax1;
  double ndt2 = cour * dx[2] / cmax2;
  double ndt3 = cour * dx[3] / cmax3;

  return (1./(1./ndt1 + 1./ndt2 + 1./ndt3));
}

/**
 * @brief Employs the flux-interpolated constrained transport scheme from Toth 2000. This smooths the fluxes with a special (diffusive) operator and maintains the no monopoles condition 
 * 
 * @param F1 Flux in direction 1 to be filled with smoothed values
 * @param F2 Flux in direction 2 to be filled with smoothed values
 * @param F3 Flux in direction 3 to be filled with smoothed values
 */
void flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3)
{
  static double emf1[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];
  static double emf2[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];
  static double emf3[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];

  // Calculate EMFs via average to corners (Toth approach)
  #pragma omp parallel
  {
    #pragma omp for collapse(3)
    // Eq. 25 in Gammie 03
    ZSLOOP(0, N1, 0, N2, 0, N3) {
      emf3[i][j][k] =  0.25*(F1[i][j][k][B2] + F1[i][j - 1][k][B2]
                           - F2[i][j][k][B1] - F2[i - 1][j][k][B1]);
      emf2[i][j][k] = -0.25*(F1[i][j][k][B3] + F1[i][j][k-1][B3]
                           - F3[i][j][k][B1] - F3[i-1][j][k][B1]);
      emf1[i][j][k] =  0.25*(F2[i][j][k][B3] + F2[i][j][k-1][B3]
                           - F3[i][j][k][B2] - F3[i][j-1][k][B2]);
    }

    // Rewrite EMFs as fluxes, after Toth
    #pragma omp for collapse(3) nowait
    ZSLOOP(0, N1, 0, N2 - 1, 0, N3 - 1) {
      F1[i][j][k][B1] = 0.;
      F1[i][j][k][B2] = 0.5 * (emf3[i][j][k] + emf3[i][j + 1][k]);
      F1[i][j][k][B3] = -0.5 * (emf2[i][j][k] + emf2[i][j][k + 1]);
    }
    #pragma omp for collapse(3) nowait
    ZSLOOP(0, N1 - 1, 0, N2, 0, N3 - 1) {
      F2[i][j][k][B1] = -0.5 * (emf3[i][j][k] + emf3[i + 1][j][k]);
      F2[i][j][k][B2] = 0.;
      F2[i][j][k][B3] = 0.5 * (emf1[i][j][k] + emf1[i][j][k + 1]);
    }
    #pragma omp for collapse(3)
    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3) {
      F3[i][j][k][B1] = 0.5 * (emf2[i][j][k] + emf2[i + 1][j][k]);
      F3[i][j][k][B2] = -0.5 * (emf1[i][j][k] + emf1[i][j + 1][k]);
      F3[i][j][k][B3] = 0.;
    }
  } // omp parallel
}

/**
 * @brief Uses the Local Lax-Friedrichs flux method to calculate the Flux at a gridpoint given the metric and the primitives at the left and right of the grid zone
 * 
 * @param P_l Primitives at the left edge
 * @param P_r Primitives at the right edge
 * @param geom ggeom struct at the center
 * @param dir direction
 * @param Flux Empty flux matrix to store results
 * @param maxSpeed Empty var to store the maximum speed
 */
void lr_to_flux(double P_l[NVAR], double P_r[NVAR], struct of_geom *geom,
  int dir, double Flux[NVAR], double *maxSpeed)
{
  struct of_state state_l, state_r;
  double F_l[NVAR], F_r[NVAR], U_l[NVAR], U_r[NVAR];
  double cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, ctop;

  if(geom->g < SMALL) {
    PLOOP Flux[ip] = 0.;
    *maxSpeed = 0.;
    return;
  }

  get_state(P_l, geom, &state_l);
  get_state(P_r, geom, &state_r);

  primtoflux(P_l, &state_l, dir, geom, F_l);
  primtoflux(P_r, &state_r, dir, geom, F_r);

  primtoflux(P_l, &state_l, 0, geom, U_l);
  primtoflux(P_r, &state_r, 0, geom, U_r);

  mhd_vchar(P_l, &state_l, geom, dir, &cmax_l, &cmin_l);
  mhd_vchar(P_r, &state_r, geom, dir, &cmax_r, &cmin_r);

  // This makes the assumption that cmax ~= cmin
  cmax = fabs(MY_MAX(MY_MAX(0., cmax_l), cmax_r));
  cmin = fabs(MY_MAX(MY_MAX(0., -cmin_l), -cmin_r));
  ctop = MY_MAX(cmax, cmin);

  PLOOP Flux[ip] = 0.5*(F_l[ip] + F_r[ip] - ctop*(U_r[ip] - U_l[ip])); // Eq. 22 of Gammie 03

  *maxSpeed = ctop;
}

