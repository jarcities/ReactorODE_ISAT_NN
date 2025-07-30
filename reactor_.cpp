#include "reactor.hpp"
#include <cmath>
#include "cantera/core.h"
#include "cantera/zerodim.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <chrono>
#include <string>
// Added includes for CVODES sensitivity analysis
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <cassert>

#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

using namespace Cantera;

namespace Gl
{

    shared_ptr<Solution> sol;    //= newSolution("h2o2.yaml", "ohmech", "none");
    shared_ptr<ThermoPhase> gas; //= sol->thermo();

    int nLayers = 5;
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

        int i1 = 0;
        int i2 = 0;

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

        //std::cerr<<A[ia[2]+3]<<std::endl;
        //std::cerr<<b[ib[1]+4]<<std::endl;
    }

}

using namespace Gl;

class ReactorODEs : public FuncEval
{
public:
    ReactorODEs(shared_ptr<Solution> sol)
    {
        m_gas = sol->thermo();
        m_kinetics = sol->kinetics();
        m_pressure = m_gas->pressure();
        m_nSpecies = m_gas->nSpecies();
        m_hbar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);
        m_nEqs = m_nSpecies + 1;
    }

    void eval(double t, double *y, double *ydot, double *p) override
    {
        double temperature = y[0];
        double *massFracs = &y[1];
        double *dTdt = &ydot[0];
        double *dYdt = &ydot[1];

        m_gas->setMassFractions_NoNorm(massFracs);
        m_gas->setState_TP(temperature, m_pressure);
        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(&m_hbar[0]);
        m_kinetics->getNetProductionRates(&m_wdot[0]);

        //energy equation
        double hdot_vol = 0;
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        *dTdt = -hdot_vol / (rho * cp);

        //species conservation equations
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }

    size_t neq() const override
    {
        return m_nEqs;
    }

    void getState(double *y) override
    {
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(&y[1]);
    }

private:
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics> m_kinetics;
    vector<double> m_hbar;
    vector<double> m_wdot;
    double m_pressure;
    size_t m_nSpecies;
    size_t m_nEqs;
};

// RHS function for CVODES sensitivity analysis
static int RHS_CVODES(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    ReactorODEs *odes = static_cast<ReactorODEs *>(user_data);
    double *y_data = NV_DATA_S(y);
    double *yd_data = NV_DATA_S(ydot);
    
    // Create temporary arrays for the eval function
    std::vector<double> y_vec(odes->neq());
    std::vector<double> ydot_vec(odes->neq());
    std::vector<double> p_vec; // dummy parameter vector
    
    // Copy data from SUNDIALS vectors
    for (size_t i = 0; i < odes->neq(); i++) {
        y_vec[i] = y_data[i];
    }
    
    // Call the existing eval function
    odes->eval((double)t, y_vec.data(), ydot_vec.data(), p_vec.data());
    
    // Copy result back to SUNDIALS vector
    for (size_t i = 0; i < odes->neq(); i++) {
        yd_data[i] = ydot_vec[i];
    }
    
    return 0;
}

void fromxhat(double x[], double ptcl[], int &nx, double rusr[])
{

    ptcl[0] = (x[0] * rusr[nx]) + rusr[0];

    for (int ii = 1; ii < nx; ii++)
    {

        ptcl[ii] = rusr[ii] * exp(-log(rusr[ii]) * x[ii]) - rusr[ii];
    }

    //for ( int ii = 0; ii < nx; ii++ ){ptcl[ii] = x[ii]*(100000.0*rusr[ii]);} //REMOVE THIS
}

void toxhat(double ptcl[], double x[], int &nx, double rusr[])
{

    x[0] = (ptcl[0] - rusr[0]) / rusr[nx];

    for (int ii = 1; ii < nx; ii++)
    {

        x[ii] = -log((ptcl[ii] + rusr[ii]) / rusr[ii]) / log(rusr[ii]);
    }

    //for ( int ii = 0; ii < nx; ii++ ){x[ii] = ptcl[ii]/(100000.0*rusr[ii]);} //REMOVE THIS
}

double fAct(double x)
{

    return x * tanh(log(1.0 + exp(x)));
}

void myfnn(int &nx, double x[], double fnn[])
{

    static int bbbb;

    double x1[100];
    double x2[100];

    if (bbbb != 7777)
    {
        Gl::initfnn();
        bbbb = 7777;
    }

    for (int ii = 0; ii < n1[0]; ii++)
    {
        x1[ii] = x[ii]; // initialize the input
    }

    for (int ll = 0; ll < nLayers; ll++)
    {

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;
            for (int jj = 0; jj < n1[ll]; jj++)
            {
                //x2[kk] += A[ ia[ll] + jj + (kk-1)*n1[ll] ]*x1[jj];
                x2[kk] += A[ia[ll] + kk + (jj - 1) * n2[ll]] * x1[jj];
            }
            x2[kk] += b[ib[ll] + kk];

            if (ll < nLayers - 1)
            {
                x2[kk] = fAct(x2[kk]);
            }
        }

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x1[kk] = x2[kk];
        }
    }

    for (int kk = 0; kk < nx; kk++)
    {
        fnn[kk] = x2[kk];
    }
}

//The actual code is put into a function that can be called from the main program
void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
           double rusr[], double f[], double g[], double h[])
{

    double Y[nx - 1];
    double T[1];
    double ptcl[nx];
    double *solution;
    double aTol = 1e-8; //rusr[2*nx];
    double rTol = 1e-8; //rusr[2*nx+1];
    double dt = rusr[2 * nx + 2];
    double dx = rusr[2 * nx + 3];
    double p = rusr[2 * nx + 4];
    double fnn[nx];

    static int aaaa;

    if (aaaa != 7777)
    {
        Gl::initfgh();
        aaaa = 7777;
    }

    fromxhat(x, ptcl, nx, rusr);

    T[0] = ptcl[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y[ii - 1] = ptcl[ii];
    }

    gas->setState_TPY(T[0], p, Y);

    //create ODE RHS evaluator
    ReactorODEs odes = ReactorODEs(sol);

    double tnow = 0.0;

    //double dt = 1e-4;

    shared_ptr<Integrator> integrator(newIntegrator("CVODE"));

    integrator->initialize(tnow, odes);

    integrator->setTolerances(aTol, rTol);

    integrator->integrate(dt);

    solution = integrator->solution();

    toxhat(solution, f, nx, rusr);

    myfnn(nx, x, fnn);

    for (int ii = 0; ii < nx; ii++)
    {
        f[ii] = f[ii] - x[ii] - fnn[ii];
    }

    //JACOBIAN START - OLD CENTRAL DIFFERENCE APPROACH (COMMENTED OUT)
    /*
    if (need[1] == 1)
    {

        double xp[nx];
        double xm[nx];
        double fp[nf];
        double fm[nf];

        for (int ii = 0; ii < nx; ii++)
        {

            for (int jj = 0; jj < nx; jj++)
            {
                xp[jj] = x[jj];
                xm[jj] = x[jj];
            }

            xp[ii] += dx;
            xm[ii] -= dx;

            tnow = 0.0;
            fromxhat(xp, ptcl, nx, rusr);
            T[0] = ptcl[0];
            for (int ii = 1; ii < nx; ii++)
            {
                Y[ii - 1] = ptcl[ii];
            }
            gas->setState_TPY(T[0], p, Y);
            integrator->initialize(tnow, odes);
            integrator->integrate(dt);
            solution = integrator->solution();
            toxhat(solution, fp, nx, rusr);
            myfnn(nx, xp, fnn);
            for (int ii = 0; ii < nx; ii++)
            {
                fp[ii] = fp[ii] - xp[ii] - fnn[ii];
            }

            tnow = 0.0;
            fromxhat(xm, ptcl, nx, rusr);
            T[0] = ptcl[0];
            for (int ii = 1; ii < nx; ii++)
            {
                Y[ii - 1] = ptcl[ii];
            }
            gas->setState_TPY(T[0], p, Y);
            integrator->initialize(tnow, odes);
            integrator->integrate(dt);
            solution = integrator->solution();
            toxhat(solution, fm, nx, rusr);
            myfnn(nx, xm, fnn);
            for (int ii = 0; ii < nx; ii++)
            {
                fm[ii] = fm[ii] - xm[ii] - fnn[ii];
            }

            for (int jj = 0; jj < nx; jj++)
            {
                g[jj + ii * (nx)] = 1.0 * (fp[jj] - fm[jj]) / (2 * dx);
            }
        }
    }
    */
    //JACOBIAN END - OLD CENTRAL DIFFERENCE APPROACH
    
    //JACOBIAN START - NEW SENSITIVITY-BASED APPROACH
    if (need[1] == 1)
    {
        // Create SUNDIALS context
        SUNContext sunctx;
        int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
        assert(flag >= 0);

        // Convert initial conditions to physical space
        fromxhat(x, ptcl, nx, rusr);
        T[0] = ptcl[0];
        for (int ii = 1; ii < nx; ii++)
        {
            Y[ii - 1] = ptcl[ii];
        }
        gas->setState_TPY(T[0], p, Y);

        // Create ReactorODEs object for sensitivity analysis
        ReactorODEs odes_sens = ReactorODEs(sol);
        size_t NEQ = odes_sens.neq();

        // Allocate state vector and set initial conditions
        N_Vector y = N_VNew_Serial(NEQ, sunctx);
        assert(y);
        double *y_data = NV_DATA_S(y);
        odes_sens.getState(y_data);

        // Create CVODES solver
        void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
        assert(cvode_mem);
        flag = CVodeInit(cvode_mem, RHS_CVODES, 0.0, y);
        assert(flag >= 0);

        // Set tolerances
        flag = CVodeSStolerances(cvode_mem, rTol, aTol);
        assert(flag >= 0);
        flag = CVodeSetUserData(cvode_mem, &odes_sens);
        assert(flag >= 0);
        flag = CVodeSetMaxNumSteps(cvode_mem, 50000);
        assert(flag >= 0);
        flag = CVodeSetMaxStep(cvode_mem, 1e-6);
        assert(flag >= 0);

        // Set linear solver
        SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
        assert(A);
        SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
        assert(LS);
        flag = CVodeSetLinearSolver(cvode_mem, LS, A);
        assert(flag >= 0);

        // Setup sensitivity analysis
        int Ns = (int)NEQ;
        N_Vector *yS = N_VCloneVectorArray(Ns, y);
        assert(yS);
        for (int is = 0; is < Ns; ++is)
        {
            for (int j = 0; j < (int)NEQ; ++j)
            {
                // Initialize sensitivity vectors as identity matrix
                NV_Ith_S(yS[is], j) = (is == j ? 1.0 : 0.0);
            }
        }

        // Initialize forward sensitivity
        flag = CVodeSensInit(cvode_mem, Ns, CV_STAGGERED, nullptr, yS);
        assert(flag >= 0);
        flag = CVodeSensEEtolerances(cvode_mem);
        assert(flag >= 0);

        // Set sensitivity parameters
        std::vector<sunrealtype> p_vec(NEQ);
        std::vector<sunrealtype> pbar_vec(NEQ);
        for (int i = 0; i < NEQ; i++)
        {
            p_vec[i] = y_data[i];
            pbar_vec[i] = (fabs(y_data[i]) > 1e-12) ? fabs(y_data[i]) : 1.0;
        }
        flag = CVodeSetSensParams(cvode_mem, p_vec.data(), pbar_vec.data(), nullptr);
        assert(flag >= 0);

        // Integrate with sensitivity
        double t_final = dt;
        double t_current = 0.0;
        flag = CVode(cvode_mem, t_final, y, &t_current, CV_NORMAL);
        if (flag < 0)
        {
            std::cerr << "CVODES integration failed with flag: " << flag << std::endl;
            // Clean up and fall back to finite difference if needed
            N_VDestroy(y);
            N_VDestroyVectorArray(yS, Ns);
            CVodeFree(&cvode_mem);
            SUNMatDestroy(A);
            SUNLinSolFree(LS);
            SUNContext_Free(&sunctx);
            return; // or implement fallback
        }

        // Get final sensitivity matrix
        flag = CVodeGetSens(cvode_mem, &t_current, yS);
        assert(flag >= 0);

        // Extract final solution in physical space
        std::vector<double> final_solution(NEQ);
        for (int i = 0; i < (int)NEQ; i++)
        {
            final_solution[i] = NV_Ith_S(y, i);
        }

        // Convert final solution to normalized space
        std::vector<double> f_sens(nx);
        toxhat(final_solution.data(), f_sens.data(), nx, rusr);

        // Apply NN correction
        myfnn(nx, x, fnn);
        for (int ii = 0; ii < nx; ii++)
        {
            f_sens[ii] = f_sens[ii] - x[ii] - fnn[ii];
        }

        // Build jacobian from sensitivity matrix
        // J[i][j] represents df_i/dx_j
        for (int j = 0; j < nx; ++j)  // columns (input variables)
        {
            double *Sj = NV_DATA_S(yS[j]);
            
            // Transform sensitivities from physical to normalized space
            std::vector<double> dSol_dxj(nx);
            toxhat(Sj, dSol_dxj.data(), nx, rusr);
            
            // Apply chain rule for transformations and NN
            for (int i = 0; i < nx; ++i)  // rows (output variables)
            {
                // The jacobian includes: d(final_state - initial_state - NN_pred)/d(initial_state)
                // = d(final_state)/d(initial_state) - I - d(NN_pred)/d(initial_state)
                double jacobian_element = dSol_dxj[i];
                
                if (i == j) {
                    jacobian_element -= 1.0; // subtract identity matrix
                }
                
                // Subtract NN jacobian contribution (finite difference for NN)
                double fnn_plus[nx], fnn_minus[nx];
                double x_plus[nx], x_minus[nx];
                for (int k = 0; k < nx; k++) {
                    x_plus[k] = x[k];
                    x_minus[k] = x[k];
                }
                x_plus[j] += dx;
                x_minus[j] -= dx;
                
                myfnn(nx, x_plus, fnn_plus);
                myfnn(nx, x_minus, fnn_minus);
                
                double dNN_dx = (fnn_plus[i] - fnn_minus[i]) / (2.0 * dx);
                jacobian_element -= dNN_dx;
                
                g[i + j * nx] = jacobian_element;
            }
        }

        // Clean up SUNDIALS objects
        N_VDestroy(y);
        N_VDestroyVectorArray(yS, Ns);
        CVodeFree(&cvode_mem);
        SUNMatDestroy(A);
        SUNLinSolFree(LS);
        SUNContext_Free(&sunctx);
    }
    //JACOBIAN END - NEW SENSITIVITY-BASED APPROACH
}

void mymix(int &nx, double x1[], double x2[], double alpha[], int iusr[], double rusr[])
{

    double Y1[nx - 1], Y2[nx - 1];
    double H1, H2;
    double T1[1], T2[1], d;
    double p = OneAtm; //rusr[2*nx+4];

    T1[0] = x1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y1[ii - 1] = x1[ii];
    }
    T2[0] = x2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y2[ii - 1] = x2[ii];
    }

    gas->setState_TPY(T1[0], p, Y1);

    H1 = gas->enthalpy_mass();

    gas->setState_TPY(T2[0], p, Y2);

    H2 = gas->enthalpy_mass();

    d = H2 - H1;
    H1 += alpha[0] * d;
    H2 -= alpha[0] * d; //mix enthalpies

    for (int ii = 0; ii < nx - 1; ii++)
    {
        d = Y2[ii] - Y1[ii];
        Y1[ii] += alpha[0] * d;
        Y2[ii] -= alpha[0] * d; //mix mass fractions
    }

    d = alpha[0] * (T2 - T1);

    gas->setState_TPY(T1[0] + d, p, Y1);
    gas->setState_HP(H1, p);

    T1[0] = gas->temperature();

    gas->setState_TPY(T2[0] - d, p, Y2);
    gas->setState_HP(H2, p);

    T2[0] = gas->temperature();

    x1[0] = T1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        x1[ii] = Y1[ii - 1];
    }
    x2[0] = T2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        x2[ii] = Y2[ii - 1];
    }
}