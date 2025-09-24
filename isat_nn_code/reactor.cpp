//#include "reactor.hpp"
//#include "custom.hpp"
//#include <cmath>
//#include <cvodes/cvodes.h>
//#include <nvector/nvector_serial.h>
//#include <sunmatrix/sunmatrix_dense.h>
//#include <sunlinsol/sunlinsol_dense.h>
//#include <sundials/sundials_types.h>
//#include <sundials/sundials_context.h>
//#include <cassert>
//#include <sundials/sundials_types.h> //ADDED

//#ifndef SUN_COMM_NULL
//#define SUN_COMM_NULL NULL
//#endif

//using namespace Cantera;

//static int RHS(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
//{
//    ReactorODEs *odes = static_cast<ReactorODEs *>(user_data);
//    double *y_data = NV_DATA_S(y);
//    double *yd_data = NV_DATA_S(ydot);
//    odes->eval((double)t, y_data, yd_data, nullptr); //nullptr for p parameter
//    return 0;
//}

//void CVODES_INTEGRATE(ReactorODEs &odes, double dt, double aTol, double rTol, double *solution)
//{
//    static SUNContext sunctx = nullptr;
//    static bool initialized = false;

//    if (!initialized)
//    {
//        int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
//        assert(flag >= 0);
//        initialized = true;
//    }

//    size_t NEQ = odes.neq();

//    N_Vector y = N_VNew_Serial(NEQ, sunctx);
//    assert(y);
//    double *y_data = NV_DATA_S(y);
//    odes.getState(y_data);

//    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
//    assert(cvode_mem);

//    double t0 = 0.0;
//    double t = t0;
//    int flag = CVodeInit(cvode_mem, RHS, t0, y);
//    assert(flag >= 0);

//    flag = CVodeSStolerances(cvode_mem, rTol, aTol);
//    assert(flag >= 0);
//    flag = CVodeSetUserData(cvode_mem, &odes);
//    assert(flag >= 0);
//    flag = CVodeSetMaxNumSteps(cvode_mem, 50000);
//    assert(flag >= 0);
//    flag = CVodeSetMaxStep(cvode_mem, 1e-6); //dt or 1e-6
//    assert(flag >= 0);

//    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
//    assert(A);
//    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
//    assert(LS);
//    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
//    assert(flag >= 0);

//    //integrate
//    flag = CVode(cvode_mem, dt, y, &t, CV_NORMAL);
//    assert(flag >= 0);

//    //copy solution
//    for (size_t i = 0; i < NEQ; i++)
//    {
//        solution[i] = y_data[i];
//    }

//    //cleanup
//    N_VDestroy(y);
//    CVodeFree(&cvode_mem);
//    SUNMatDestroy(A);
//    SUNLinSolFree(LS);
//}

///////////////////////////////////////////////////////////////////////////
//void CVODES_INTEGRATE_WITH_SENS(ReactorODEs &odes, double dt, double aTol, double rTol, double *solution, bool SENS)
//{
//    static SUNContext sunctx = nullptr;
//    static bool initialized = false;

//    if (!initialized)
//    {
//        int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
//        assert(flag >= 0);
//        initialized = true;
//    }

//    size_t NEQ = odes.neq();

//    N_Vector y = N_VNew_Serial(NEQ, sunctx);
//    assert(y);
//    double *y_data = NV_DATA_S(y);
//    odes.getState(y_data);

//    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
//    assert(cvode_mem);

//    double t0 = 0.0;
//    double t = t0;
//    int flag = CVodeInit(cvode_mem, RHS, t0, y);
//    assert(flag >= 0);

//    flag = CVodeSStolerances(cvode_mem, rTol, aTol);
//    assert(flag >= 0);
//    flag = CVodeSetUserData(cvode_mem, &odes);
//    assert(flag >= 0);
//    flag = CVodeSetMaxNumSteps(cvode_mem, 50000);
//    assert(flag >= 0);
//    flag = CVodeSetMaxStep(cvode_mem, 1e-6); //dt or 1e-6
//    assert(flag >= 0);

//    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
//    assert(A);
//    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
//    assert(LS);
//    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
//    assert(flag >= 0);

//    //set sensitivity stuff
//    if (SENS)
//    {
//        size_t NS = odes.neq();
//        size_t NEQ = odes.neq();
//        y = N_VNew_Serial(NEQ, sunctx);
//        N_Vector* yS = N_VCloneVectorArray(NS, y);
//        assert(yS);
//        sunrealtype *ySdata = N_VGetArrayPointer(yS);
//        flag = CVodeSensInit(cvode_mem, NS, CV_SIMULTANEOUS, NULL, yS);
//        assert(flag >= 0);
//        flag = CVodeSetSensErrCon(cvode_mem, SUNTRUE);
//        assert(flag >= 0);
//        flag = CVodeSensEEtolerances(cvode_mem);
//        assert(flag >= 0);
//        flag = CVodeSetSensParams(cvode_mem, NULL, NULL, NULL);
//        assert(flag >= 0);
//    }

//    //integrate
//    flag = CVode(cvode_mem, dt, y, &t, CV_NORMAL);
//    assert(flag >= 0);

//    //get sensitivity
//    if (SENS)
//    {
//        flag = CVodeGetSens(cvode_mem, &t, yS);
//    }

//    //copy solution
//    for (size_t i = 0; i < NEQ; i++)
//    {
//        solution[i] = y_data[i];
//    }

//    //cleanup
//    N_VDestroy(y);
//    if (SENS)
//    {
//        N_VDestroyVectorArray(yS, NS);
//    }
//    CVodeFree(&cvode_mem);
//    SUNMatDestroy(A);
//    SUNLinSolFree(LS);
//}
///////////////////////////////////////////////////////////////////////////

//namespace Gl
//{

//    shared_ptr<Solution> sol;    //= newSolution("h2o2.yaml", "ohmech", "none");
//    shared_ptr<ThermoPhase> gas; //= sol->thermo();

//    int nLayers = 6;   //number of MLP layers
//    int nNeurons = 10; //number of neurons in each hidden layer
//    int nx = 11;       //number of input/output variables

//    int ia[100];
//    int ib[100];
//    int n1[100];
//    int n2[100];

//    double A[1000000];
//    double b[10000]; //lines 15 to 21 are work variables for reading in the MLP weights

//    void initfgh()
//    {

//        sol = newSolution("h2o2.yaml", "ohmech", "none"); //initialize the Cantera gas object,
//        //the inputs can be modified to change the chemical mechanism
//        gas = sol->thermo();
//    }

//    void initfnn()
//    { //initialize the weights for the f^{MLP} function

//        int i1 = 0;
//        int i2 = 0; //work variables for reading in the MLP weights

//        char file1[50];
//        char file2[50];

//        for (int ll = 1; ll <= nLayers; ll++)
//        {

//            ia[ll - 1] = i1;
//            ib[ll - 1] = i2;

//            sprintf(file1, "./A%d.csv", ll);
//            sprintf(file2, "./B%d.csv", ll);

//            n1[ll - 1] = nNeurons;
//            n2[ll - 1] = nNeurons;
//            if (ll == 1)
//            {
//                n1[ll - 1] = nx;
//            }
//            if (ll == nLayers)
//            {
//                n2[ll - 1] = nx;
//            }

//            FILE *pFile;
//            float a;
//            pFile = fopen(file1, "r+");
//            for (int ii = 0; ii < n1[ll - 1] * n2[ll - 1]; ii++)
//            {
//                fscanf(pFile, "%f", &a);
//                A[i1] = a;
//                i1 = i1 + 1;
//            }
//            fclose(pFile);

//            pFile = fopen(file2, "r+");
//            for (int ii = 0; ii < n2[ll - 1]; ii++)
//            {
//                fscanf(pFile, "%f", &a);
//                b[i2] = a;
//                i2 = i2 + 1;
//            }
//            fclose(pFile);
//        }
//    }

//}

//using namespace Gl;

////"custom.hpp" is the ReactorODEs class from lines 21-147 of the file "custom.cpp" which can be found
////at this URL: https://cantera.org/3.1/examples/cxx/custom.html (retrieved 06/05/2025)

//void fromxhat(double x[], double ptcl[], int &nx, double rusr[])
//{
//    //this function converts the normalized vector x into temperature and mass fractions for one particle
//    //x[] is the input, ptcl[] is the output, nx indicates the number of dimensions of both x and ptcl
//    //rusr[] are user-supplied normalization variables

//    ptcl[0] = (x[0] * rusr[nx]) + rusr[0]; //ptcl[0] is the temperature, in K

//    for (int ii = 1; ii < nx; ii++)
//    {

//        ptcl[ii] = rusr[ii] * exp(-log(rusr[ii]) * x[ii]) - rusr[ii]; //ptcl[ii] is the mass fraction of the ii-th
//                                                                      //species in the chemical mechanism
//    }
//}

//void toxhat(double ptcl[], double x[], int &nx, double rusr[])
//{
//    //this function converts a particle's temperature and mass fractions into the normalized vector x
//    //x[] is the input, ptcl[] is the output, nx indicates the number of dimensions of both x and ptcl
//    //rusr[] are user-supplied normalization variables

//    x[0] = (ptcl[0] - rusr[0]) / rusr[nx]; //x[0] is the normalized temperature

//    for (int ii = 1; ii < nx; ii++)
//    {

//        x[ii] = -log((ptcl[ii] + rusr[ii]) / rusr[ii]) / log(rusr[ii]); //x[ii] is the normalized mass fraction
//                                                                        //of the ii-th species in the chemical mechanism
//    }
//}

//double fAct(double x)
//{ //activation function of the hidden layers,
//    //here a Mish function is used

//    return x * tanh(log(1.0 + exp(x)));
//}

//void myfnn(int &nx, double x[], double fnn[])
//{
//    //this function evaluates f^{}

//    static int bbbb; //dummy variable used to call "initfnn" the first time "myfnn" is called

//    double x1[100];
//    double x2[100]; //work arrays

//    if (bbbb != 7777)
//    {
//        Gl::initfnn();
//        bbbb = 7777;
//    } //if "myfnn" is called for the first time, initialize the
//      //f^{MLP} data structure by reading in the weights

//    for (int ii = 0; ii < n1[0]; ii++)
//    {
//        x1[ii] = x[ii]; //initialize the input
//    }

//    for (int ll = 0; ll < nLayers; ll++)
//    {

//        for (int kk = 0; kk < n2[ll]; kk++)
//        {
//            x2[kk] = 0.0;
//            for (int jj = 0; jj < n1[ll]; jj++)
//            {
//                x2[kk] += A[ia[ll] + jj + kk * n1[ll]] * x1[jj]; //apply weights in a dense layer
//            }
//            x2[kk] += b[ib[ll] + kk]; //apply the bias in a dense layer

//            if (ll < nLayers - 1)
//            {
//                x2[kk] = fAct(x2[kk]); //apply the activation function in the hidden layers
//            }
//        }

//        for (int kk = 0; kk < n2[ll]; kk++)
//        {
//            x1[kk] = x2[kk];
//        }
//    }

//    for (int kk = 0; kk < nx; kk++)
//    {
//        fnn[kk] = 1.0 * (x2[kk]);
//    } //pass the output
//}

////myfgh is the function passed to ISAT
//void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
//           double rusr[], double f[], double g[], double h[])
//{

//    double Y[nx - 1];             //mass fraction
//    double T[1];                  //temperature
//    double ptcl[nx];              //particle properties
//    double aTol = 1e-8;           //rusr[2*nx];
//    double rTol = 1e-8;           //rusr[2*nx+1]; //absolute and relative tolerances for the ODE integrator
//    double dt = rusr[2 * nx + 2]; //time step over which to integrate
//    double dx = rusr[2 * nx + 3]; //spatial increment in x for Jacobian evaluation
//    double p = rusr[2 * nx + 4];  //user-specified pressure
//    int mode = iusr[0];
//    double fnn[nx]; //f^{MLP}

//    static int aaaa; //if "myfgh" is called for the first time, initialize the
//                     //myfgh data structure by creating the Cantera solution object

//    if (aaaa != 7777)
//    {
//        Gl::initfgh();
//        aaaa = 7777;
//    } //initialize "myfgh" on the first call

//    fromxhat(x, ptcl, nx, rusr); //convert from normalized x to particle properties

//    T[0] = ptcl[0];
//    for (int ii = 1; ii < nx; ii++)
//    {
//        Y[ii - 1] = ptcl[ii];
//    } //extract temperature and mass fractions

//    gas->setState_TPY(T[0], p, Y); //set the gas state, using the particle's temperature and mass fractions,
//    //and the user-specified pressure

//    ReactorODEs odes = ReactorODEs(sol); //create the ODE RHS evaluator

//    double tnow = 0.0; //initialize time

//    /////////////////////////////////////////////////////////
//    double solution_arr[nx];
//    CVODES_INTEGRATE(odes, dt, aTol, rTol, solution_arr);
//    /////////////////////////////////////////////////////////

//    toxhat(solution_arr, f, nx, rusr); //normalize the gas properties

//    if (mode == 2)
//    {
//        myfnn(nx, x, fnn); //evaluate f^{MLP}
//        for (int ii = 0; ii < nx; ii++)
//        {
//            f[ii] = f[ii] - x[ii] - fnn[ii];
//        }
//    }
//    //f(x) is the increment of x minus f^{MLP}(x)
//    else
//    {
//        for (int ii = 0; ii < nx; ii++)
//        {
//            f[ii] = f[ii] - x[ii];
//        }
//    }

//    if (need[1] == 1)
//    {
//        //TODO: implement CVODES_INTEGRATE_WITH_SENS to get g mode = 3, mode != 2
//    }
//    //if (need[1] == 1)
//    //{ //this block is called when a Jacobian is needed

//    //    double xp[nx];
//    //    double xm[nx];
//    //    double fp[nf];
//    //    double fm[nf];

//    //    for (int ii = 0; ii < nx; ii++)
//    //    {

//    //        for (int jj = 0; jj < nx; jj++)
//    //        {
//    //            xp[jj] = x[jj];
//    //            xm[jj] = x[jj];
//    //        }

//    //        xp[ii] += dx;
//    //        xm[ii] -= dx;

//    //        tnow = 0.0;
//    //        fromxhat(xp, ptcl, nx, rusr);
//    //        T[0] = ptcl[0];
//    //        for (int ii = 1; ii < nx; ii++)
//    //        {
//    //            Y[ii - 1] = ptcl[ii];
//    //        }
//    //        gas->setState_TPY(T[0], p, Y);
//    //        ReactorODEs odes_p = ReactorODEs(sol);
//    //        /////////////////////////////////////////////////////////
//    //        double solution_arr_p[nx];
//    //        CVODES_INTEGRATE(odes_p, dt, aTol, rTol, solution_arr_p);
//    //        /////////////////////////////////////////////////////////
//    //        toxhat(solution_arr_p, fp, nx, rusr);
//    //        if (mode == 2)
//    //        {
//    //            myfnn(nx, xp, fnn);
//    //            for (int ii = 0; ii < nx; ii++)
//    //            {
//    //                fp[ii] = fp[ii] - xp[ii] - fnn[ii];
//    //            }
//    //        }
//    //        //evaluate f(x) at x+dx in the jj+1-th entry of x
//    //        else
//    //        {
//    //            for (int ii = 0; ii < nx; ii++)
//    //            {
//    //                fp[ii] = fp[ii] - xp[ii];
//    //            }
//    //        }

//    //        tnow = 0.0;
//    //        fromxhat(xm, ptcl, nx, rusr);
//    //        T[0] = ptcl[0];
//    //        for (int ii = 1; ii < nx; ii++)
//    //        {
//    //            Y[ii - 1] = ptcl[ii];
//    //        }
//    //        gas->setState_TPY(T[0], p, Y);
//    //        ReactorODEs odes_m = ReactorODEs(sol);
//    //        /////////////////////////////////////////////////////////
//    //        double solution_arr_m[nx];
//    //        CVODES_INTEGRATE(odes_m, dt, aTol, rTol, solution_arr_m);
//    //        /////////////////////////////////////////////////////////
//    //        toxhat(solution_arr_m, fm, nx, rusr);
//    //        if (mode == 2)
//    //        {
//    //            myfnn(nx, xm, fnn);
//    //            for (int ii = 0; ii < nx; ii++)
//    //            {
//    //                fm[ii] = fm[ii] - xm[ii] - fnn[ii];
//    //            }
//    //        }
//    //        //evaluate f(x) at x-dx in the jj+1-th entry of x
//    //        else
//    //        {
//    //            for (int ii = 0; ii < nx; ii++)
//    //            {
//    //                fm[ii] = fm[ii] - xm[ii];
//    //            }
//    //        }

//    //        for (int jj = 0; jj < nx; jj++)
//    //        {
//    //            g[jj + ii * (nx)] = 1.0 * (fp[jj] - fm[jj]) / (2 * dx);
//    //        }
//    //        //calculate the jj+1-th partial derivative
//    //    }
//    //}
//}

//void mymix(int &nx, double ptcl1[], double ptcl2[], double alpha[], int iusr[], double rusr[])
//{
//    //mix two particles, conserving mass and energy

//    double Y1[nx - 1], Y2[nx - 1]; //mass fractions
//    double H1, H2;                 //enthalpies
//    double T1[1], T2[1];           //temperatures
//    double d;                      //work variable
//    double p = OneAtm;             //rusr[2*nx+4];

//    T1[0] = ptcl1[0];
//    for (int ii = 1; ii < nx; ii++)
//    {
//        Y1[ii - 1] = ptcl1[ii];
//    }
//    T2[0] = ptcl2[0];
//    for (int ii = 1; ii < nx; ii++)
//    {
//        Y2[ii - 1] = ptcl2[ii];
//    }
//    //extract temperature and mass fractios for the two particles

//    gas->setState_TPY(T1[0], p, Y1); //initialize the gas to the state of the first particle

//    H1 = gas->enthalpy_mass(); //get the enthalpy of the first particle

//    gas->setState_TPY(T2[0], p, Y2); //initialize the gas to the state of the second particle

//    H2 = gas->enthalpy_mass(); //get the enthalpy of the second particle

//    d = H2 - H1;
//    H1 += alpha[0] * d;
//    H2 -= alpha[0] * d; //mix enthalpies
//    //alpha is amount by which to mix the particles (0 is no mixing, 0.5 is complete mixing)

//    for (int ii = 0; ii < nx - 1; ii++)
//    {
//        d = Y2[ii] - Y1[ii];
//        Y1[ii] += alpha[0] * d;
//        Y2[ii] -= alpha[0] * d; //mix mass fractions
//    }

//    d = alpha[0] * (T2 - T1);

//    gas->setState_TPY(T1[0] + d, p, Y1);
//    gas->setState_HP(H1, p);

//    T1[0] = gas->temperature();

//    gas->setState_TPY(T2[0] - d, p, Y2);
//    gas->setState_HP(H2, p);

//    T2[0] = gas->temperature(); //set the particle's thermodynamic states to the new mixed values, and
//    //extract the corresponding temperatures

//    ptcl1[0] = T1[0];
//    for (int ii = 1; ii < nx; ii++)
//    {
//        ptcl1[ii] = Y1[ii - 1];
//    }
//    ptcl2[0] = T2[0];
//    for (int ii = 1; ii < nx; ii++)
//    {
//        ptcl2[ii] = Y2[ii - 1];
//    } //pass the new temperatures and mass fractions to the output
//}

#include "reactor.hpp"
#include "custom.hpp"
#include <cassert>
#include <cmath>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <vector>

#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

using namespace Cantera;

//////////////////////
struct CVUserData
{
    ReactorODEs *odes;
    sunrealtype *p;
    int NP;
};
//////////////////////

///////////////////////////////////////////////////////////////////////////
static int RHS(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    auto *ud = static_cast<CVUserData *>(user_data);
    ud->odes->eval((double)t, NV_DATA_S(y), NV_DATA_S(ydot), ud->p);
    return 0;
}
///////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
void CVODES_INTEGRATE(ReactorODEs &odes, double dt, double aTol, double rTol,
                      double *solution, double *JAC /* = nullptr */)
{
    static SUNContext sunctx = nullptr;
    static bool initialized = false;
    if (!initialized)
    {
        int f = SUNContext_Create(SUN_COMM_NULL, &sunctx);
        assert(f >= 0);
        initialized = true;
    }

    //num of equations and solution vector
    const sunindextype NEQ = (sunindextype)odes.neq();
    N_Vector y = N_VNew_Serial(NEQ, sunctx);
    assert(y);
    odes.getState(NV_DATA_S(y));

    //create/init SUNDIAL solver and set tolerances
    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
    assert(cvode_mem);
    int flag = CVodeInit(cvode_mem, RHS, 0.0, y);
    assert(flag >= 0);
    flag = CVodeSStolerances(cvode_mem, rTol, aTol);
    assert(flag >= 0);

    //set user data (ode, params, num params)
    CVUserData ud{&odes, nullptr, 0};
    flag = CVodeSetUserData(cvode_mem, &ud);
    assert(flag >= 0);

    //set max time steps
    flag = CVodeSetMaxNumSteps(cvode_mem, 50000);
    assert(flag >= 0);
    flag = CVodeSetMaxStep(cvode_mem, dt);
    assert(flag >= 0);

    //set linear solver (dense data structures)
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    assert(A);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    assert(LS);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    assert(flag >= 0);

    //init sens stuff
    N_Vector *yS = nullptr;
    int NS = 0;

    //dummy variables
    std::vector<sunrealtype> pvec;
    std::vector<sunrealtype> pbar;

    if (JAC)
    {
        NS = (int)NEQ;

        //1) seed sensitivity ICs as identity (IC sensitivities)
        yS = N_VCloneVectorArray(NS, y);
        assert(yS);
        for (int i = 0; i < NS; ++i)
        {
            N_VConst(0.0, yS[i]);
            NV_Ith_S(yS[i], i) = 1.0;
        }

        flag = CVodeSensInit(cvode_mem, NS, CV_SIMULTANEOUS, /*fS*/ NULL, yS);
        assert(flag >= 0);

        //2) provide a dummy parameter vector so internal DQ path is valid
        pvec.assign(NS, 0.0); //values don't matter (RHS ignores p)
        pbar.assign(NS, 1.0); //scale (positive)
        ud.p = pvec.data();
        ud.NP = NS;
        flag = CVodeSetUserData(cvode_mem, &ud);
        assert(flag >= 0);
        flag = CVodeSetSensParams(cvode_mem, ud.p, pbar.data(), /*plist*/ nullptr);
        assert(flag >= 0);

        //3) enable sensitivities with internal DQ (no custom fS)

        flag = CVodeSensEEtolerances(cvode_mem);
        assert(flag >= 0);
        flag = CVodeSetSensErrCon(cvode_mem, SUNTRUE);
        assert(flag >= 0);
    }

    //integrate once
    double t = 0.0;
    flag = CVode(cvode_mem, dt, y, &t, CV_NORMAL);
    assert(flag >= 0);

    //final state
    double *y_data = NV_DATA_S(y);
    for (sunindextype i = 0; i < NEQ; ++i)
        solution[i] = y_data[i];

    //copy sensitivities (column-major): JAC[j + i*NEQ] = ∂y_j(t)/∂y0_i
    if (JAC && NS > 0)
    {
        flag = CVodeGetSens(cvode_mem, &t, yS);
        assert(flag >= 0);
        for (int i = 0; i < NS; ++i)
            for (sunindextype j = 0; j < NEQ; ++j)
                JAC[j + i * NEQ] = NV_Ith_S(yS[i], j);
    }

    //cleanup
    if (yS)
        N_VDestroyVectorArray(yS, NS);
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
}
////////////////////////////////////////////////////////////////////////////////////////

namespace Gl
{
    shared_ptr<Solution> sol;
    shared_ptr<ThermoPhase> gas;

    int nLayers = 6;
    int nNeurons = 10;
    int nx = 11;

    int ia[100];
    int ib[100];
    int n1[100];
    int n2[100];

    double A[1000000];
    double b[10000];

    void initfgh()
    {
        sol = newSolution("h2o2.yaml", "ohmech", "none");
        gas = sol->thermo();
    }

    void initfnn()
    {
        int i1 = 0, i2 = 0;
        char file1[50], file2[50];

        for (int ll = 1; ll <= nLayers; ll++)
        {
            ia[ll - 1] = i1;
            ib[ll - 1] = i2;

            sprintf(file1, "./A%d.csv", ll);
            sprintf(file2, "./B%d.csv", ll);

            n1[ll - 1] = (ll == 1) ? nx : nNeurons;
            n2[ll - 1] = (ll == nLayers) ? nx : nNeurons;

            FILE *pFile;
            float a;

            pFile = fopen(file1, "r+");
            for (int ii = 0; ii < n1[ll - 1] * n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                A[i1++] = a;
            }
            fclose(pFile);

            pFile = fopen(file2, "r+");
            for (int ii = 0; ii < n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                b[i2++] = a;
            }
            fclose(pFile);
        }
    }
}

using namespace Gl;

//helpers (same formulas you already had)
void fromxhat(double x[], double ptcl[], int &nx, double rusr[])
{
    ptcl[0] = (x[0] * rusr[nx]) + rusr[0];
    for (int ii = 1; ii < nx; ii++)
    {
        ptcl[ii] = rusr[ii] * std::exp(-std::log(rusr[ii]) * x[ii]) - rusr[ii];
    }
}

void toxhat(double ptcl[], double x[], int &nx, double rusr[])
{
    x[0] = (ptcl[0] - rusr[0]) / rusr[nx];
    for (int ii = 1; ii < nx; ii++)
    {
        x[ii] = -std::log((ptcl[ii] + rusr[ii]) / rusr[ii]) / std::log(rusr[ii]);
    }
}

double fAct(double x) { return x * std::tanh(std::log(1.0 + std::exp(x))); }

void myfnn(int &nx, double x[], double fnn[])
{
    static int bbbb;
    double x1[100], x2[100];

    if (bbbb != 7777)
    {
        Gl::initfnn();
        bbbb = 7777;
    }

    for (int ii = 0; ii < n1[0]; ii++)
        x1[ii] = x[ii];

    for (int ll = 0; ll < nLayers; ll++)
    {
        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;
            for (int jj = 0; jj < n1[ll]; jj++)
            {
                x2[kk] += A[ia[ll] + jj + kk * n1[ll]] * x1[jj];
            }
            x2[kk] += b[ib[ll] + kk];
            if (ll < nLayers - 1)
                x2[kk] = fAct(x2[kk]);
        }
        for (int kk = 0; kk < n2[ll]; kk++)
            x1[kk] = x2[kk];
    }
    for (int kk = 0; kk < nx; kk++)
        fnn[kk] = x2[kk];
}

//myfgh is the function passed to ISAT
void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
           double rusr[], double f[], double g[], double h[])
{
    //double Y[nx - 1];             //mass fraction
    //double T[1];                  //temperature
    //double ptcl[nx];              //particle properties
    double aTol = 1e-8;           //rusr[2*nx];
    double rTol = 1e-8;           //rusr[2*nx+1]; //absolute and relative tolerances for the ODE integrator
    double dt = rusr[2 * nx + 2]; //time step over which to integrate
    //double dx = rusr[2 * nx + 3]; //spatial increment in x for Jacobian evaluation
    //double p = rusr[2 * nx + 4];  //user-specified pressure
    int mode = iusr[0];
    //double fnn[nx]; //f^{MLP}
    std::vector<double> solution_arr(nx, 0.0);
    std::vector<double> JAC; 

    //init first
    static int aaaa;
    if (aaaa != 7777)
    {
        Gl::initfgh();
        aaaa = 7777;
    }

    //ic -> normalize
    std::vector<double> ptcl(nx);
    fromxhat(x, ptcl.data(), nx, rusr);

    //get T and Y
    double T = ptcl[0];
    std::vector<double> Y(nx - 1);
    for (int i = 1; i < nx; ++i)
        Y[i - 1] = ptcl[i];

    //set state and build RHS
    double p = rusr[2 * nx + 4];
    gas->setState_TPY(T, p, Y.data());
    ReactorODEs odes = ReactorODEs(sol);

    ////integrate ONCE; ask for sensitivities only if Jacobian is needed
    //double aTol = 1e-8;           //or rusr[2*nx]
    //double rTol = 1e-8;           //or rusr[2*nx+1]
    //double dt = rusr[2 * nx + 2]; //integration horizon
    //int mode = iusr[0];           //your mode flag
    
    double *JAC_ptr = nullptr;
    if (need[1] == 1)
    {
        JAC.resize(nx * nx, 0.0);
        JAC_ptr = JAC.data();
    }

    //integrate
    CVODES_INTEGRATE(odes, dt, aTol, rTol, solution_arr.data(), JAC_ptr);
    toxhat(solution_arr.data(), f, nx, rusr);

    if (mode == 2)
    {
        double fnn[100];
        myfnn(nx, x, fnn);
        for (int i = 0; i < nx; ++i)
            f[i] = f[i] - x[i] - fnn[i];
    }
    else
    {
        for (int i = 0; i < nx; ++i)
            f[i] = f[i] - x[i];
    }

    //If Jacobian requested: g = B * JAC * A - I
    if (need[1] == 1)
    {
        //A = d(y0)/d(x) at ICs (ptcl from fromxhat)
        std::vector<double> A_diag(nx);
        A_diag[0] = rusr[nx]; //dT/dx0
        for (int i = 1; i < nx; ++i)
            A_diag[i] = -(ptcl[i] + rusr[i]) * std::log(rusr[i]); //dYi/dx_i

        //B = d(x)/d(y(t)) at final state (solution_arr)
        std::vector<double> B_diag(nx);
        B_diag[0] = 1.0 / rusr[nx]; //dx0/dT
        for (int i = 1; i < nx; ++i)
            B_diag[i] = -1.0 / ((solution_arr[i] + rusr[i]) * std::log(rusr[i])); //dx_i/dYi

        //assemble (column-major)
        for (int col = 0; col < nx; ++col)
        {
            const double Acol = A_diag[col];
            for (int row = 0; row < nx; ++row)
            {
                double val = B_diag[row] * JAC[row + col * nx] * Acol;
                if (row == col)
                    val -= 1.0; //-I from f(...)-x
                g[row + col * nx] = val;
            }
        }
    }
}

void mymix(int &nx, double ptcl1[], double ptcl2[], double alpha[], int iusr[], double rusr[])
{
    //conserve mass and energy
    std::vector<double> Y1(nx - 1), Y2(nx - 1);
    double T1 = ptcl1[0], T2 = ptcl2[0];
    for (int i = 1; i < nx; ++i)
    {
        Y1[i - 1] = ptcl1[i];
        Y2[i - 1] = ptcl2[i];
    }

    double p = OneAtm; //or rusr[2*nx+4]
    gas->setState_TPY(T1, p, Y1.data());
    double H1 = gas->enthalpy_mass();
    gas->setState_TPY(T2, p, Y2.data());
    double H2 = gas->enthalpy_mass();

    double dH = H2 - H1;
    H1 += alpha[0] * dH;
    H2 -= alpha[0] * dH;

    for (int i = 0; i < nx - 1; ++i)
    {
        double dY = Y2[i] - Y1[i];
        Y1[i] += alpha[0] * dY;
        Y2[i] -= alpha[0] * dY;
    }

    double dT = alpha[0] * (T2 - T1);

    gas->setState_TPY(T1 + dT, p, Y1.data());
    gas->setState_HP(H1, p);
    T1 = gas->temperature();

    gas->setState_TPY(T2 - dT, p, Y2.data());
    gas->setState_HP(H2, p);
    T2 = gas->temperature();

    ptcl1[0] = T1;
    for (int i = 1; i < nx; ++i)
        ptcl1[i] = Y1[i - 1];
    ptcl2[0] = T2;
    for (int i = 1; i < nx; ++i)
        ptcl2[i] = Y2[i - 1];
}