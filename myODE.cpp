#include <cvodes/cvodes.h> /* prototypes for CVODES fcts., consts. */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include "canteraFxns.hpp"

#define Ith(v, i) NV_Ith_S(v, i - 1) /* i-th vector component i=1..NEQ */
#define IJth(A, i, j) \
  SM_ELEMENT_D(A, i - 1, j - 1) /* (i,j)-th matrix component i,j=1..NEQ */



/* Problem Constants */

#define NEQ   nSp+1               /* number of equations  */
#define Y1    SUN_RCONST(1.0) /* initial y components */
#define Y2    SUN_RCONST(0.0)
#define Y3    SUN_RCONST(0.0)
#define RTOL  SUN_RCONST(1.0e-9) /* scalar relative tolerance            */
#define ATOL1 SUN_RCONST(1.0e-8) /* vector absolute tolerance components */
#define ATOL2 SUN_RCONST(1.0e-12)
#define ATOL3 SUN_RCONST(1.0e-8)
#define T0    SUN_RCONST(0.0)  /* initial time           */
#define T1    SUN_RCONST(0.000004)  /* first output time      */
#define TMULT SUN_RCONST(5.0) /* output time factor     */
#define NOUT  1               /* number of output times */

#define NS nSp+1 /* number of sensitivities computed */

#define ZERO SUN_RCONST(0.0)
#define ONE SUN_RCONST(1.0)

/* Type : UserData */

using namespace Cantera;

typedef struct
{
  sunrealtype pressure;  // pressure of the constant pressure reactor
  sunrealtype p[NEQ]; /* problem parameters */
  shared_ptr< ThermoPhase > gas;
  shared_ptr< Kinetics > kinetics;
}* UserData;

/* Functions Called by the Solver */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fS(int Ns, sunrealtype t, N_Vector y, N_Vector ydot, int iS,
              N_Vector yS, N_Vector ySdot, void* user_data, N_Vector tmp1,
              N_Vector tmp2);

static int ewt(N_Vector y, N_Vector w, void* user_data);

/* Private functions to output results */

static void PrintOutput(void* cvode_mem, sunrealtype t, N_Vector u);
static void PrintOutputS(N_Vector* uS);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char* argv[], sunbooleantype* sensi,
                        int* sensi_meth, sunbooleantype* err_con);

static void WrongArgs(char* name);

/* Private function to check function return values */

static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char* argv[])
{
	
  double MW[nSp];
  double hbar[nSp];
  initCantera2(MW, hbar);
	
  SUNContext sunctx;
  sunrealtype t, tout;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, iout;
  UserData data;
  FILE* FID;
  char fname[256];

  sunrealtype pbar[NS];
  int is;
  N_Vector* yS;
  sunbooleantype sensi, err_con;
  int sensi_meth;

  data      = NULL;
  y         = NULL;
  yS        = NULL;
  A         = NULL;
  LS        = NULL;
  cvode_mem = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* User data structure */
  data = (UserData)malloc(sizeof *data);
  if (check_retval((void*)data, "malloc", 2)) { return (1); }

  /* Initialize sensitivity variables (reaction rates for this problem) */
  for (int ii=0; ii<NS; ii++){data->p[ii]=SUN_RCONST(1.0);}
  data->pressure = SUN_RCONST(101325.0);
  //data->gas = Gl::gas;
  //data->kinetics = Gl::kinetics;

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* Initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

  sunrealtype p[NEQ];
  for (int ii=0; ii<NEQ; ii++){p[ii]=data->p[ii];}

  /* Initialize y */
  for ( int ii = 0; ii < NS; ii++ ){
	  Ith(y,ii+1) = ZERO;
  }  
  
  Ith(y,1) = SUN_RCONST(1500.0);
  Ith(y,2) = SUN_RCONST(0.099999);
  Ith(y,3) = SUN_RCONST(0.000001);
  Ith(y,5) = SUN_RCONST(0.9);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  /* Call CVodeWFtolerances to specify a user-supplied function ewt that sets
   * the multiplicative error weights w_i for use in the weighted RMS norm */
  retval = CVodeWFtolerances(cvode_mem, ewt);
  if (check_retval(&retval, "CVodeWFtolerances", 1)) { return (1); }

  /* Attach user data */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, NULL);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }

  printf(" \n3-species kinetics problem\n");

  /* Sensitivity-related settings */
  if (sensi)
  {
	  
	for ( int ii = 0; ii < NS; ii++ ){
		pbar[ii] = 0.000001;
	}
	
	pbar[0] = 1.0;
    
    /* Set sensitivity initial conditions */
    yS = N_VCloneVectorArray(NS, y);
    if (check_retval((void*)yS, "N_VCloneVectorArray", 0)) { return (1); }
    for (is = 0; is < NS; is++) { N_VConst(ZERO, yS[is]); }
	
	for ( int ii = 0; ii < NS; ii++ ){
		Ith(yS[ii],ii+1) = ONE;
	}
	
	
	//sunrealtype *ySdata = N_VGetArrayPointer(yS);	
	//printf("%12.4e %12.4e %12.4e \n", ySdata[0], ySdata[1], ySdata[2]);

    /* Call CVodeSensInit1 to activate forward sensitivity computations
     * and allocate internal memory for COVEDS related to sensitivity
     * calculations. Computes the right-hand sides of the sensitivity
     * ODE, one at a time */
    retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, yS);
    if (check_retval(&retval, "CVodeSensInit", 1)) { return (1); }

    /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
     * variables based on the rolerances supplied for states variables and
     * the scaling factor pbar */
    retval = CVodeSensEEtolerances(cvode_mem);
    if (check_retval(&retval, "CVodeSensEEtolerances", 1)) { return (1); }

    /* Set sensitivity analysis optional inputs */
    /* Call CVodeSetSensErrCon to specify the error control strategy for
     * sensitivity variables */
    retval = CVodeSetSensErrCon(cvode_mem, err_con);
    if (check_retval(&retval, "CVodeSetSensErrCon", 1)) { return (1); }

    /* Call CVodeSetSensParams to specify problem parameter information for
     * sensitivity calculations */
    retval = CVodeSetSensParams(cvode_mem, data->p, pbar, NULL);
    if (check_retval(&retval, "CVodeSetSensParams", 1)) { return (1); }

    printf("Sensitivity: YES ");
    if (sensi_meth == CV_SIMULTANEOUS) { printf("( SIMULTANEOUS +"); }
    else if (sensi_meth == CV_STAGGERED) { printf("( STAGGERED +"); }
    else { printf("( STAGGERED1 +"); }
    if (err_con) { printf(" FULL ERROR CONTROL )"); }
    else { printf(" PARTIAL ERROR CONTROL )"); }
  }
  else { printf("Sensitivity: NO "); }

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */

  printf("\n\n");
  printf("===========================================");
  printf("============================\n");
  printf("     T     Q       H      NST           y1");
  printf("           y2           y3    \n");
  printf("===========================================");
  printf("============================\n");

  retval = CVodeSetMaxNumSteps(cvode_mem, 5000000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) { return (1); }

  for (iout = 1, tout = T1; iout <= NOUT; iout++, tout *= TMULT)
  {
	retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) { break; }
    PrintOutput(cvode_mem, t, y);

    /* Call CVodeGetSens to get the sensitivity solution vector after a
     * successful return from CVode */
    if (sensi)
    {
      retval = CVodeGetSens(cvode_mem, &t, yS);
      if (check_retval(&retval, "CVodeGetSens", 1)) { break; }
      PrintOutputS(yS);
    }
    printf("-----------------------------------------");
    printf("------------------------------\n");
  }

  /* Print final statistics to the screen */
  printf("\nFinal Statistics:\n");
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* Print final statistics to a file in CSV format */
  strcpy(fname, "cvsRoberts_FSA_dns_stats");
  if (sensi)
  {
    if (sensi_meth == CV_SIMULTANEOUS) { strcat(fname, "_-sensi_sim"); }
    else if (sensi_meth == CV_STAGGERED) { strcat(fname, "_-sensi_stg"); }
    else { strcat(fname, "_-sensi_stg1"); }
    if (err_con) { strcat(fname, "_t"); }
    else { strcat(fname, "_f"); }
  }
  strcat(fname, ".csv");
  FID    = fopen(fname, "w");
  retval = CVodePrintAllStats(cvode_mem, FID, SUN_OUTPUTFORMAT_CSV);
  fclose(FID);

  /* Free memory */
  N_VDestroy(y); /* Free y vector */
  if (sensi) { N_VDestroyVectorArray(yS, NS); /* Free yS vector */ }
  free(data);               /* Free user data */
  CVodeFree(&cvode_mem);    /* Free CVODES memory */
  SUNLinSolFree(LS);        /* Free the linear solver memory */
  SUNMatDestroy(A);         /* Free the matrix memory */
  SUNContext_Free(&sunctx); /* Free the SUNDIALS context */

  return (0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype Y[NEQ];
  UserData data;
  sunrealtype pressure,hdot,wdot[nSp],rho,cp,hbar[nSp];
  
  double sumddy, ddy[nSp];

  data = (UserData)user_data;
  pressure   = Gl::gas->pressure();
  
  int nsp = Gl::gas->nSpecies();

  for (int ii=0; ii < NEQ; ii++){Y[ii] = Ith(y,ii+1);}  
  
  
  Gl::gas->setMassFractions(&Y[1]);
  Gl::gas->setState_TP(Y[0],pressure);
  Gl::gas->getPartialMolarEnthalpies(&hbar[0]);
  Gl::kinetics->getNetProductionRates(&wdot[0]);
  
  hdot = 0.0;
  for (int ii = 0; ii < nSp; ii++){hdot+=hbar[ii]*wdot[ii];}
  
  rho = Gl::gas->density();
  cp = Gl::gas->cp_mass();
  
  Ith(ydot,1) = -hdot/(rho*cp);
  
  sumddy = 0.0;
  
  for (int ii = 0; ii < nSp; ii++){
	  Ith(ydot,ii+2) = (wdot[ii]*Gl::gas->molecularWeight(ii))/rho;
	  ddy[ii] = (wdot[ii]*Gl::gas->molecularWeight(ii))/rho;
	  sumddy += ddy[ii];
  }
  
  return (0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /*sunrealtype dy;
  double fyp[NEQ], fym[NEQ], yp[NEQ], ym[NEQ];
  UserData data;
  sunrealtype p1, p2, p3;
  int rett;
  
  fyp = N_VNew_Serial(NEQ, sunctx);
  fym = N_VNew_Serial(NEQ, sunctx);
  yp = N_VNew_Serial(NEQ, sunctx);
  ym = N_VNew_Serial(NEQ, sunctx);
  
  dy = 0.000000001;
  
  for (int jj = 0; jj < NEQ; jj++){
	  
	  for (int ii = 0; ii < NEQ; ii++){
		Ith(yp,ii+1) = Ith(y,ii+1); Ith(ym,ii+1) = Ith(y,ii+1);  
	  }
	  
	  Ith(yp, jj+1) = Ith(yp, jj+1) + dy; 
	  Ith(ym, jj+1) = Ith(ym, jj+1) - dy;
	  
	  rett = f(t,yp,fyp, user_data);
	  rett = f(t,ym,fym, user_data);
	  
	  for (int ii=0; ii<NEQ; ii++){
		 IJth(J, ii+1, jj+1) = (Ith(fyp,ii+1) - Ith(fym,ii+1))/(2*dy); 
	  }	  
  }*/

  return (0);
}

/*
 * fS routine. Compute sensitivity r.h.s.
 */

static int fS(int Ns, sunrealtype t, N_Vector y, N_Vector ydot, int iS,
              N_Vector yS, N_Vector ySdot, void* user_data, N_Vector tmp1,
              N_Vector tmp2)
{
  SUNMatrix myJ;
  
  int rett;
  sunrealtype ysdot1;
  
  rett = Jac(t, y, ydot, myJ, user_data, ydot, ydot, ydot);
  
  for (int ii = 0; ii < NS; ii++){
	  ysdot1 = ZERO;
	  for (int jj = 0; jj < NS; jj++){
		ysdot1 += IJth( myJ, ii+1, jj+1)*Ith(yS, jj+1); 
	  }
	  Ith(ySdot, ii+1) = ysdot1;
  }
    
  return (0);
}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void* user_data)
{
  int i;
  sunrealtype yy, ww, rtol, atol[NEQ];

  rtol    = RTOL;
  atol[0] = ATOL1;  
  for ( int ii=1; ii<NEQ; ii++ ){atol[ii] = ATOL2;}
  
  for (i = 1; i <= NEQ; i++)
  {
    yy = Ith(y, i);
    ww = rtol * abs(yy) + atol[i-1];
    if (ww <= 0.0) { return (-1); }
    Ith(w, i) = 1.0 / ww;
  }

  return (0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */

static void ProcessArgs(int argc, char* argv[], sunbooleantype* sensi,
                        int* sensi_meth, sunbooleantype* err_con)
{
  *sensi      = SUNFALSE;
  *sensi_meth = -1;
  *err_con    = SUNFALSE;

  if (argc < 2) { WrongArgs(argv[0]); }

  if (strcmp(argv[1], "-nosensi") == 0) { *sensi = SUNFALSE; }
  else if (strcmp(argv[1], "-sensi") == 0) { *sensi = SUNTRUE; }
  else { WrongArgs(argv[0]); }

  if (*sensi)
  {
    if (argc != 4) { WrongArgs(argv[0]); }

    if (strcmp(argv[2], "sim") == 0) { *sensi_meth = CV_SIMULTANEOUS; }
    else if (strcmp(argv[2], "stg") == 0) { *sensi_meth = CV_STAGGERED; }
    else if (strcmp(argv[2], "stg1") == 0) { *sensi_meth = CV_STAGGERED1; }
    else { WrongArgs(argv[0]); }

    if (strcmp(argv[3], "t") == 0) { *err_con = SUNTRUE; }
    else if (strcmp(argv[3], "f") == 0) { *err_con = SUNFALSE; }
    else { WrongArgs(argv[0]); }
  }
}

static void WrongArgs(char* name)
{
  printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n", name);
  printf("         sensi_meth = sim, stg, or stg1\n");
  printf("         err_con    = t or f\n");

  exit(0);
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(void* cvode_mem, sunrealtype t, N_Vector u)
{
  long int nst;
  int qu, retval;
  sunrealtype hu, *udata;

  udata = N_VGetArrayPointer(u);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

  printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%14.6Le %14.6Le %14.6Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%14.6e %14.6e %14.6e \n", udata[0], udata[1], udata[2]);
#else
  printf("%14.6e %14.6e %14.6e \n", udata[0], udata[1], udata[2]);
#endif
}

/*
 * Print sensitivities.
*/

static void PrintOutputS(N_Vector* uS)
{
  sunrealtype* sdata;

  sdata = N_VGetArrayPointer(uS[0]);
  printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[1]);
  printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[2]);
  printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
