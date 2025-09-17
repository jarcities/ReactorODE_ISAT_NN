#include "reactor_.hpp"
#include "custom_.hpp"
#include <cmath>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes.h>  // Use cvodes instead of cvode for sensitivity analysis
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cassert>
#include <iostream>
#include <vector>

#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif


using namespace Cantera;

namespace Gl
{

    shared_ptr<Solution> sol;    // = newSolution("h2o2.yaml", "ohmech", "none");
    shared_ptr<ThermoPhase> gas; // = sol->thermo();

    int nLayers = 6;   // number of MLP layers
    int nNeurons = 10; // number of neurons in each hidden layer
    int nx = 11;       // number of input/output variables

    int ia[100];
    int ib[100];
    int n1[100];
    int n2[100];

    double A[1000000];
    double b[10000]; // lines 15 to 21 are work variables for reading in the MLP weights

    void initfgh()
    {

        sol = newSolution("h2o2.yaml", "ohmech", "none"); // initialize the Cantera gas object,
        // the inputs can be modified to change the chemical mechanism
        gas = sol->thermo();
    }

    void initfnn()
    { // initialize the weights for the f^{MLP} function

        int i1 = 0;
        int i2 = 0; // work variables for reading in the MLP weights

        char file1[50];
        char file2[50];

        for (int ll = 1; ll <= nLayers; ll++)
        {

            ia[ll - 1] = i1;
            ib[ll - 1] = i2;

            sprintf(file1, "./A%d.csv", ll);
            sprintf(file2, "./B%d.csv", ll);

            n1[ll - 1] = nNeurons;
            n2[ll - 1] = nNeurons;
            if (ll == 1)
            {
                n1[ll - 1] = nx;
            }
            if (ll == nLayers)
            {
                n2[ll - 1] = nx;
            }

            FILE *pFile;
            float a;
            pFile = fopen(file1, "r+");
            for (int ii = 0; ii < n1[ll - 1] * n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                A[i1] = a;
                i1 = i1 + 1;
            }
            fclose(pFile);

            pFile = fopen(file2, "r+");
            for (int ii = 0; ii < n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                b[i2] = a;
                i2 = i2 + 1;
            }
            fclose(pFile);
        }
    }

}

using namespace Gl;


// "custom.hpp" is the ReactorODEs class from lines 21-147 of the file "custom.cpp" which can be found
// at this URL: https://cantera.org/3.1/examples/cxx/custom.html (retrieved 06/05/2025)

void fromxhat(double x[], double ptcl[], int &nx, double rusr[])
{
    ptcl[0] = (x[0] * rusr[nx]) + rusr[0];

    for (int ii = 1; ii < nx; ii++)
    {

        ptcl[ii] = rusr[ii] * exp(-log(rusr[ii]) * x[ii]) - rusr[ii]; 
    }
}

void toxhat(double ptcl[], double x[], int &nx, double rusr[])
{
    x[0] = (ptcl[0] - rusr[0]) / rusr[nx]; 

    for (int ii = 1; ii < nx; ii++)
    {

        x[ii] = -log((ptcl[ii] + rusr[ii]) / rusr[ii]) / log(rusr[ii]);
    }
}

double fAct(double x)
{ // activation function of the hidden layers,
    // here a Mish function is used

    return x * tanh(log(1.0 + exp(x)));
}

////////////////////////////////////////////////////////////////////////////////////////
double dfAct(double x)
{
    // Mish(x) = x * tanh(softplus(x)),  softplus = log1p(exp(x))
    double s = log1p(exp(x));
    double t = tanh(s);
    double sech2 = 1.0 - t*t;
    double sig = 1.0 / (1.0 + exp(-x)); // d softplus / dx
    return t + x * sech2 * sig;
}

////////////////////////////////////////////////////////////////////////////////////////

void myfnn(int need[], int &nx, double x[], double fnn[], double jnn[])
{
    // this function evaluates f^{}

    static int bbbb; // dummy variable used to call "initfnn" the first time "myfnn" is called

    double x1[100];
    double x2[100];         // work arrays

    ////////////////////////////////////////////
    double jx[100]; 
    double jac[100][100]; 
    double jtemp[100][100]; 
    ////////////////////////////////////////////

    if (bbbb != 7777)
    {
        Gl::initfnn();
        bbbb = 7777;
    } 
      

    ////////////////////////////////////////////
    if (need[1] == 1)
    {
        for (int i = 0; i < nx; ++i) 
        { 
            for (int j = 0; j < nx; ++j) 
            { 
                jtemp[i][j] = (i == j) ? 1.0 : 0.0; 
            } 
        } 
    }
    ////////////////////////////////////////////

    for (int ii = 0; ii < n1[0]; ii++)
    {
        x1[ii] = x[ii]; // initialize the input
    }

    for (int ll = 0; ll < nLayers; ll++)
    {

        // layer propagation
        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;
            for (int jj = 0; jj < n1[ll]; jj++)
            {
                x2[kk] += A[ia[ll] + jj + kk * n1[ll]] * x1[jj]; // apply weights in a dense layer
            }
            x2[kk] += b[ib[ll] + kk]; // apply the bias in a dense layer

            if (ll < nLayers - 1)
            {
                ///////////////////////
                if (need[1] == 1)
                {
                    jx[kk] = dfAct(x2[kk]); 
                }   
                ///////////////////////
                x2[kk] = fAct(x2[kk]);  // apply the activation function in the hidden layers
            }
        }

        ////////////////////////////////////////////////////////////
        if (need[1] == 1)
        {
            for (int i = 0; i < n2[ll]; ++i) 
            { 
                for (int j = 0; j < nx; ++j) 
                { 
                    double sum = 0.0; 
                    for (int k = 0; k < n1[ll]; ++k) 
                        sum += A[ia[ll] + k + i * n1[ll]] * jtemp[k][j]; 
                    jac[i][j] = (ll < nLayers - 1 ? jx[i] * sum : sum); 
                    jtemp[i][j] = jac[i][j]; 
                } 
            } 
        }
        ////////////////////////////////////////////////////////////

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x1[kk] = x2[kk];
        }
    }

    for (int kk = 0; kk < nx; kk++)
    {
        fnn[kk] = 1.0 * (x2[kk]);

        //////////////////////////////////////
        if (need[1] == 1)
        {
            for (int jj = 0; jj < nx; ++jj) 
            {
                jnn[jj * nx + kk] = jtemp[kk][jj];
            }
        }
        //////////////////////////////////////
    } 
}

static int RHS(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    ReactorODEs *odes = static_cast<ReactorODEs *>(user_data);
    double *y_data = NV_DATA_S(y); 
    double *yd_data = NV_DATA_S(ydot);
    odes->eval((double)t, y_data, yd_data, nullptr);
    return 0; 
} 

void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
           double rusr[], double f[], double g[], double h[])
{

    double Y[nx - 1];             // mass fraction
    double T[1];                  // temperature
    double ptcl[nx];              // particle properties
    double *solution;             // Cantera solution object
    double aTol = 1e-4;           // rusr[2*nx];
    double rTol = 1e-4;           // rusr[2*nx+1]; //absolute and relative tolerances for the ODE integrator
    double dt = rusr[2 * nx + 2]; // time step over which to integrate
    double dx = rusr[2 * nx + 3]; // spatial increment in x for Jacobian evaluation
    double p = rusr[2 * nx + 4];  // user-specified pressure
    int mode = iusr[0];
    double fnn[nx]; 
    static std::vector<double> jnn(nx * nx, 0.0); 

    //sensitivity stuff
    double tnow = 0.0;
    double t = tnow;
    const double tfinal = dt;

    static int aaaa;
    int flag;
    static SUNContext sunctx; 
    static bool cvodes_init = false;
    
    // Static variables for CVODES reuse - MAJOR PERFORMANCE IMPROVEMENT
    static void *m_cvode_mem = nullptr;
    static N_Vector y_static = nullptr;
    static N_Vector *yS_static = nullptr;
    static SUNMatrix A_static = nullptr;
    static SUNLinearSolver LS_static = nullptr;
    static int Ns_static = 0;
    static bool cvodes_objects_init = false; 

    if (aaaa != 7777)
    {
        Gl::initfgh();
        aaaa = 7777;
    } 

    if (!cvodes_init)
    {
        flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
        assert(flag == 0);
        cvodes_init = true;
    }

    fromxhat(x, ptcl, nx, rusr); 

    T[0] = ptcl[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y[ii - 1] = ptcl[ii];
    } 

    gas->setState_TPY(T[0], p, Y);
    ReactorODEs odes = ReactorODEs(sol);
    size_t NEQ = odes.neq();
    assert((int)NEQ == nx && "Mismatch between reactor state size and normalized vector size");

    // Initialize CVODES objects only once - HUGE PERFORMANCE GAIN
    if (!cvodes_objects_init) {
        y_static = N_VNew_Serial(NEQ, sunctx);
        assert(y_static);
        m_cvode_mem = CVodeCreate(CV_BDF, sunctx);
        assert(m_cvode_mem);
        A_static = SUNDenseMatrix(NEQ, NEQ, sunctx);
        assert(A_static);
        LS_static = SUNLinSol_Dense(y_static, A_static, sunctx);
        assert(LS_static);
        
        // Initial setup of CVODES
        double *y_temp = NV_DATA_S(y_static);
        odes.getState(y_temp);
        flag = CVodeInit(m_cvode_mem, RHS, tnow, y_static);
        assert(flag >= 0);
        flag = CVodeSetLinearSolver(m_cvode_mem, LS_static, A_static);
        assert(flag >= 0);
        
        // Always setup sensitivity arrays for maximum flexibility
        Ns_static = (int)NEQ;
        yS_static = N_VCloneVectorArray(Ns_static, y_static);
        assert(yS_static);
        // Initialize sensitivity vectors
        for (int is = 0; is < Ns_static; ++is)
        {
            for (int j = 0; j < (int)NEQ; ++j)
            {
                NV_Ith_S(yS_static[is], j) = (is == j ? 1.0 : 0.0);
            }
        }
        flag = CVodeSensInit(m_cvode_mem, Ns_static, CV_SIMULTANEOUS, nullptr, yS_static); //CV_STAGGERED or CV_SIMULTANEOUS
        assert(flag >= 0);
        
        cvodes_objects_init = true;
    }

    // Reset and reinitialize for this call
    double *y_data = NV_DATA_S(y_static);
    odes.getState(y_data);
    
    // Reinitialize with new initial conditions (much faster than full setup)
    flag = CVodeReInit(m_cvode_mem, tnow, y_static);
    assert(flag >= 0);
    flag = CVodeSStolerances(m_cvode_mem, rTol, aTol);
    assert(flag >= 0);
    CVodeSetSensDQMethod(m_cvode_mem, CV_CENTERED, 1e-4); //CV_CENTERED or CV_FORWARD
    flag = CVodeSetUserData(m_cvode_mem, &odes);
    assert(flag >= 0);
    flag = CVodeSetMaxNumSteps(m_cvode_mem, 50000);
    assert(flag >= 0);
    flag = CVodeSetMaxStep(m_cvode_mem, dt); //1e-6 original
    assert(flag >= 0);

    //jacobian
    if (need[1] == 1)
    {
        // Reinitialize sensitivity vectors
        for (int is = 0; is < Ns_static; ++is)
        {
            for (int j = 0; j < (int)NEQ; ++j)
            {
                NV_Ith_S(yS_static[is], j) = (is == j ? 1.0 : 0.0);
            }
        }
        //forward sens - reinitialize with existing vectors
        flag = CVodeSensReInit(m_cvode_mem, CV_SIMULTANEOUS, yS_static); //CV_STAGGERED or CV_SIMULTANEOUS
        assert(flag >= 0);
        flag = CVodeSetSensErrCon(m_cvode_mem, SUNFALSE); //SUNTRUE or SUNFALSE (sensitivity does not control integrator)
        assert(flag >= 0);
        sunrealtype s_rTol = 1e-5;                
        std::vector<sunrealtype> s_aTol(Ns_static, 1e-5); 
        flag = CVodeSensSStolerances(m_cvode_mem, s_rTol, s_aTol.data());
        // flag = CVodeSensEEtolerances(m_cvode_mem);
        assert(flag >= 0);
        std::vector<sunrealtype> p(Ns_static, 0.0);
        std::vector<sunrealtype> pbar(Ns_static, 1.0);
        flag = CVodeSetSensParams(m_cvode_mem, p.data(), pbar.data(), nullptr);
        assert(flag >= 0);
    } else {
        // Turn off sensitivity if not needed for this call
        flag = CVodeSensToggleOff(m_cvode_mem);
        assert(flag >= 0);
    }

    //integrate
    flag = CVode(m_cvode_mem, tfinal, y_static, &t, CV_NORMAL);
    if (flag < 0)
    {
        std::cout << "integration failed -> " << flag << std::endl;
        return;  // No cleanup needed - static objects persist
    }

    toxhat(NV_DATA_S(y_static), f, nx, rusr); 

    if (mode == 2)
    {
        myfnn(need, nx, x, fnn, jnn.data()); 

        for (int ii = 0; ii < nx; ii++)
        {
            f[ii] = f[ii] - x[ii] - fnn[ii];
        }
    }
    else
    {
        for (int ii = 0; ii < nx; ii++)
        {
            f[ii] = f[ii] - x[ii];
        }
    }

    //jacobian
    if (need[1] == 1)
    {
        flag = CVodeGetSens(m_cvode_mem, &t, yS_static);
        assert(flag >= 0);
        std::vector<double> J_to_tf(nx), J_from_0(nx);

        double* y_final = NV_DATA_S(y_static);
        J_to_tf[0] = 1.0 / rusr[nx]; 
        for (int i = 1; i < nx; ++i) {
            double ri = rusr[i];
            J_to_tf[i] = -1.0 / (std::log(ri) * (y_final[i] + ri));
        }

        J_from_0[0] = rusr[nx]; 
        for (int j = 1; j < nx; ++j) {
            double ri = rusr[j];
            double y0_plus = ri * std::exp(-std::log(ri) * x[j]);
            J_from_0[j] = - y0_plus * std::log(ri);
        }

        for (int j = 0; j < (int)NEQ; ++j) {           
            double* Sj = NV_DATA_S(yS_static[j]);             
            for (int i = 0; i < (int)NEQ; ++i) {      
                double Sxhat_ij = J_to_tf[i] * Sj[i] * J_from_0[j];
                double nn_ij = (mode == 2) ? jnn[i + j*nx] : 0.0;  
                g[i + j*nx] = Sxhat_ij - (i == j ? 1.0 : 0.0) - nn_ij;
            }
        }
    }

}

void mymix(int &nx, double ptcl1[], double ptcl2[], double alpha[], int iusr[], double rusr[])
{

    double Y1[nx - 1], Y2[nx - 1]; // mass fractions
    double H1, H2;                 // enthalpies
    double T1[1], T2[1];           // temperatures
    double d;                      // work variable
    double p = OneAtm;             // rusr[2*nx+4];

    T1[0] = ptcl1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y1[ii - 1] = ptcl1[ii];
    }
    T2[0] = ptcl2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y2[ii - 1] = ptcl2[ii];
    }

    gas->setState_TPY(T1[0], p, Y1); // initialize the gas to the state of the first particle

    H1 = gas->enthalpy_mass(); // get the enthalpy of the first particle

    gas->setState_TPY(T2[0], p, Y2); // initialize the gas to the state of the second particle

    H2 = gas->enthalpy_mass(); // get the enthalpy of the second particle

    d = H2 - H1;
    H1 += alpha[0] * d;
    H2 -= alpha[0] * d; // mix enthalpies

    for (int ii = 0; ii < nx - 1; ii++)
    {
        d = Y2[ii] - Y1[ii];
        Y1[ii] += alpha[0] * d;
        Y2[ii] -= alpha[0] * d; // mix mass fractions
    }

    d = alpha[0] * (T2[0] - T1[0]);

    gas->setState_TPY(T1[0] + d, p, Y1);
    gas->setState_HP(H1, p);

    T1[0] = gas->temperature();

    gas->setState_TPY(T2[0] - d, p, Y2);
    gas->setState_HP(H2, p);

    T2[0] = gas->temperature(); // set the particle's thermodynamic states to the new mixed values, and

    ptcl1[0] = T1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        ptcl1[ii] = Y1[ii - 1];
    }
    ptcl2[0] = T2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        ptcl2[ii] = Y2[ii - 1];
    }
}
